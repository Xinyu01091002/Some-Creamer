#include <chrono>
#include <algorithm>
#include <complex>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>

#ifdef CREAMER_HAVE_OPENMP
#include <omp.h>
#endif

namespace fs = std::filesystem;
using cdouble = std::complex<double>;

template <typename T>
std::vector<T> read_vector(const fs::path& path, std::size_t count) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("Could not open input: " + path.string());
    std::vector<T> data(count);
    in.read(reinterpret_cast<char*>(data.data()), static_cast<std::streamsize>(count * sizeof(T)));
    if (!in) throw std::runtime_error("Could not read expected bytes: " + path.string());
    return data;
}

std::vector<cdouble> read_complex_vector(const fs::path& path, std::size_t count) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("Could not open input: " + path.string());
    std::vector<cdouble> data(count);
    for (std::size_t i = 0; i < count; ++i) {
        double re = 0.0, im = 0.0;
        in.read(reinterpret_cast<char*>(&re), sizeof(double));
        in.read(reinterpret_cast<char*>(&im), sizeof(double));
        if (!in) throw std::runtime_error("Could not read complex data: " + path.string());
        data[i] = cdouble(re, im);
    }
    return data;
}

template <typename T>
void write_vector(const fs::path& path, const std::vector<T>& data) {
    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("Could not open output: " + path.string());
    out.write(reinterpret_cast<const char*>(data.data()), static_cast<std::streamsize>(data.size() * sizeof(T)));
}

void write_complex_vector(const fs::path& path, const std::vector<cdouble>& data) {
    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("Could not open output: " + path.string());
    for (const auto& z : data) {
        const double re = z.real();
        const double im = z.imag();
        out.write(reinterpret_cast<const char*>(&re), sizeof(double));
        out.write(reinterpret_cast<const char*>(&im), sizeof(double));
    }
}

struct Plan {
    std::int64_t ny = 0;
    std::int64_t nx = 0;
    std::int64_t n_total = 0;
    std::int64_t n_active = 0;
    std::int64_t n_steps = 0;
    bool preserve_mean = true;
    double dx = 1.0;
    double dy = 1.0;
    double g = 9.81;
    double energy_fraction = 0.99;
    double min_active_modes = 0.0;
    double max_active_modes = 0.0;
    std::vector<std::int64_t> active;
    std::vector<std::int64_t> dest;
    std::vector<double> bz;
    std::vector<double> bphi;
    std::vector<double> d;
    std::vector<std::int64_t> order;
};

double kernel_d(double theta1, double theta2, double theta3, double g) {
    const double den_base = theta1 * (theta2 + theta3 - theta1)
                          + theta2 * (theta1 + theta3 - theta2)
                          + theta3 * (theta1 + theta2 - theta3);
    if (std::abs(den_base) < 1e-14 || theta1 == 0.0 || theta2 == 0.0 || theta3 == 0.0) {
        return 0.0;
    }
    const double num = (theta1 + theta2 + theta3)
                     * (theta1 + theta2 - theta3)
                     * (theta1 - theta2 - theta3)
                     * (theta1 + theta3 - theta2);
    return num / (12.0 * g * den_base);
}

double kernel_b(double theta1, double theta2, double theta3) {
    const double den_base = theta1 * (theta2 + theta3 - theta1)
                          + theta2 * (theta1 + theta3 - theta2)
                          + theta3 * (theta1 + theta2 - theta3);
    if (std::abs(den_base) < 1e-14) return 0.0;
    const double num = (theta1 - theta2) * (theta1 - theta2) * (theta1 + theta2)
                     - theta3 * theta3 * (2.0 * theta3 - theta1 - theta2);
    return num / (2.0 * den_base);
}

std::int64_t fft_mode_number(std::int64_t idx, std::int64_t n) {
    if (n == 1) return 0;
    if (idx <= n / 2) return idx;
    return idx - n;
}

std::int64_t positive_mod(std::int64_t a, std::int64_t n) {
    const std::int64_t r = a % n;
    return r < 0 ? r + n : r;
}

std::int64_t mode_to_index(std::int64_t m, std::int64_t n) {
    std::int64_t m_wrapped = positive_mod(m + n / 2, n) - n / 2;
    std::int64_t idx = m_wrapped;
    if (idx < 0) idx += n;
    return idx;
}

std::vector<double> build_kmag(std::int64_t ny, std::int64_t nx, double dy, double dx) {
    const double lx = static_cast<double>(nx) * dx;
    const double ly = static_cast<double>(ny) * dy;
    std::vector<double> kmag(static_cast<std::size_t>(ny * nx), 0.0);
    for (std::int64_t col = 0; col < nx; ++col) {
        const double kx = (nx == 1) ? 0.0 : (2.0 * M_PI / lx) * static_cast<double>(fft_mode_number(col, nx));
        for (std::int64_t row = 0; row < ny; ++row) {
            const double ky = (ny == 1) ? 0.0 : (2.0 * M_PI / ly) * static_cast<double>(fft_mode_number(row, ny));
            kmag[static_cast<std::size_t>(row + col * ny)] = std::hypot(kx, ky);
        }
    }
    return kmag;
}

void build_active_modes(Plan& p, const std::vector<cdouble>& zeta) {
    std::vector<std::pair<double, std::int64_t>> energy;
    energy.reserve(static_cast<std::size_t>(p.n_total));
    double total_energy = 0.0;
    for (std::int64_t idx = 0; idx < p.n_total; ++idx) {
        const double e = std::norm(zeta[static_cast<std::size_t>(idx)]);
        energy.emplace_back(e, idx);
        total_energy += e;
    }
    std::stable_sort(energy.begin(), energy.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.first > rhs.first;
    });

    std::int64_t keep_count = p.n_total;
    if (total_energy > 0.0) {
        double cumulative = 0.0;
        for (std::int64_t i = 0; i < p.n_total; ++i) {
            cumulative += energy[static_cast<std::size_t>(i)].first;
            if (cumulative / total_energy >= p.energy_fraction) {
                keep_count = i + 1;
                break;
            }
        }
    }
    if (std::isfinite(p.min_active_modes)) {
        keep_count = std::max(keep_count, std::min(static_cast<std::int64_t>(p.min_active_modes), p.n_total));
    }
    if (std::isfinite(p.max_active_modes)) {
        keep_count = std::min(keep_count, std::max<std::int64_t>(0, static_cast<std::int64_t>(p.max_active_modes)));
    }
    keep_count = std::clamp<std::int64_t>(keep_count, 0, p.n_total);

    std::vector<unsigned char> active_mask(static_cast<std::size_t>(p.n_total), 0);
    for (std::int64_t i = 0; i < keep_count; ++i) {
        active_mask[static_cast<std::size_t>(energy[static_cast<std::size_t>(i)].second)] = 1;
    }
    if (p.preserve_mean && std::abs(zeta[0]) > 0.0) {
        active_mask[0] = 1;
    } else {
        active_mask[0] = 0;
    }

    p.active.clear();
    for (std::int64_t idx = 0; idx < p.n_total; ++idx) {
        if (active_mask[static_cast<std::size_t>(idx)]) p.active.push_back(idx);
    }
    p.n_active = static_cast<std::int64_t>(p.active.size());
}

void build_interaction_plan(Plan& p, const std::vector<double>& kmag) {
    const std::size_t na = static_cast<std::size_t>(p.n_active);
    const std::size_t npair = na * na;
    p.dest.assign(npair, 0);
    p.bz.assign(npair, 0.0);
    p.bphi.assign(npair, 0.0);
    p.d.assign(npair, 0.0);

    std::vector<std::int64_t> mx(na), my(na);
    std::vector<double> active_kmag(na, 0.0);
    for (std::size_t a = 0; a < na; ++a) {
        const std::int64_t idx = p.active[a];
        const std::int64_t row = idx % p.ny;
        const std::int64_t col = idx / p.ny;
        mx[a] = fft_mode_number(col, p.nx);
        my[a] = fft_mode_number(row, p.ny);
        active_kmag[a] = kmag[static_cast<std::size_t>(idx)];
    }

#pragma omp parallel for schedule(static)
    for (std::int64_t b64 = 0; b64 < p.n_active; ++b64) {
        const std::size_t b = static_cast<std::size_t>(b64);
        for (std::size_t a = 0; a < na; ++a) {
            const std::size_t off = a + b * na;
            const std::int64_t col_k = mode_to_index(mx[a] + mx[b], p.nx);
            const std::int64_t row_k = mode_to_index(my[a] + my[b], p.ny);
            const std::int64_t dest = row_k + col_k * p.ny;
            const double kmag_k = kmag[static_cast<std::size_t>(dest)];
            p.dest[off] = dest;
            p.d[off] = kernel_d(kmag_k, active_kmag[a], active_kmag[b], p.g);
            p.bz[off] = kernel_b(active_kmag[a], active_kmag[b], kmag_k);
            p.bphi[off] = kernel_b(kmag_k, active_kmag[a], active_kmag[b]);
        }
    }

    p.order.resize(npair);
    for (std::size_t i = 0; i < npair; ++i) p.order[i] = static_cast<std::int64_t>(i);
    std::stable_sort(p.order.begin(), p.order.end(), [&](std::int64_t lhs, std::int64_t rhs) {
        return p.dest[static_cast<std::size_t>(lhs)] < p.dest[static_cast<std::size_t>(rhs)];
    });
}

Plan read_job_metadata(const fs::path& dir) {
    auto meta = read_vector<std::int64_t>(dir / "meta.bin", 5);
    auto params = read_vector<double>(dir / "params.bin", 6);
    Plan p;
    p.ny = meta[0];
    p.nx = meta[1];
    p.n_total = meta[2];
    p.n_steps = meta[3];
    p.preserve_mean = (meta[4] != 0);
    p.dx = params[0];
    p.dy = params[1];
    p.g = params[2];
    p.energy_fraction = params[3];
    p.min_active_modes = params[4];
    p.max_active_modes = params[5];

    if (p.n_total <= 0 || p.n_steps <= 0 || p.nx <= 0 || p.ny <= 0 || p.dx <= 0.0 || p.dy <= 0.0) {
        throw std::runtime_error("Invalid metadata in job folder");
    }
    return p;
}

void rhs(const Plan& p,
         const std::vector<cdouble>& zeta,
         const std::vector<cdouble>& phi,
         std::vector<cdouble>& dz,
         std::vector<cdouble>& dp) {
    std::fill(dz.begin(), dz.end(), cdouble(0.0, 0.0));
    std::fill(dp.begin(), dp.end(), cdouble(0.0, 0.0));

    const std::size_t na = static_cast<std::size_t>(p.n_active);
    const std::size_t npair = na * na;
    const std::size_t block_size = 16384;
    const std::size_t n_blocks = (npair + block_size - 1) / block_size;
    std::vector<std::vector<std::int64_t>> block_dest(n_blocks);
    std::vector<std::vector<cdouble>> block_z(n_blocks);
    std::vector<std::vector<cdouble>> block_p(n_blocks);

#pragma omp parallel for schedule(dynamic)
    for (std::int64_t block = 0; block < static_cast<std::int64_t>(n_blocks); ++block) {
        const std::size_t start = static_cast<std::size_t>(block) * block_size;
        const std::size_t stop = std::min(npair, start + block_size);
        auto& out_dest = block_dest[static_cast<std::size_t>(block)];
        auto& out_z = block_z[static_cast<std::size_t>(block)];
        auto& out_p = block_p[static_cast<std::size_t>(block)];
        out_dest.reserve(block_size);
        out_z.reserve(block_size);
        out_p.reserve(block_size);

        std::size_t cursor = start;
        while (cursor < stop) {
            const std::int64_t off0 = p.order[cursor];
            const std::int64_t dest = p.dest[static_cast<std::size_t>(off0)];
            cdouble sum_z(0.0, 0.0);
            cdouble sum_p(0.0, 0.0);

            do {
                const std::size_t off = static_cast<std::size_t>(p.order[cursor]);
                const std::size_t a = off % na;
                const std::size_t b = off / na;
                const auto za = zeta[static_cast<std::size_t>(p.active[a])];
                const auto pa = phi[static_cast<std::size_t>(p.active[a])];
                const auto zb = zeta[static_cast<std::size_t>(p.active[b])];
                const auto pb = phi[static_cast<std::size_t>(p.active[b])];
                sum_z += 3.0 * p.d[off] * pa * pb + p.bz[off] * za * zb;
                sum_p += -2.0 * p.bphi[off] * za * pb;
                ++cursor;
            } while (cursor < stop && p.dest[static_cast<std::size_t>(p.order[cursor])] == dest);

            out_dest.push_back(dest);
            out_z.push_back(sum_z);
            out_p.push_back(sum_p);
        }
    }

    for (std::size_t block = 0; block < n_blocks; ++block) {
        for (std::size_t i = 0; i < block_dest[block].size(); ++i) {
            const std::size_t k = static_cast<std::size_t>(block_dest[block][i]);
            dz[k] += block_z[block][i];
            dp[k] += block_p[block][i];
        }
    }
}

void axpy(std::vector<cdouble>& out,
          const std::vector<cdouble>& base,
          double h,
          const std::vector<cdouble>& k) {
    const std::size_t n = out.size();
    for (std::size_t i = 0; i < n; ++i) out[i] = base[i] + h * k[i];
}

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            std::cerr << "Usage: creamer_flow <job_dir>\n";
            return 2;
        }

        const fs::path dir(argv[1]);
        const auto t0 = std::chrono::steady_clock::now();
        Plan p = read_job_metadata(dir);
        const std::size_t n = static_cast<std::size_t>(p.n_total);
        const double h = -1.0 / static_cast<double>(p.n_steps);

        auto zeta = read_complex_vector(dir / "zeta0.bin", n);
        auto phi = read_complex_vector(dir / "phi0.bin", n);
        const auto kmag = build_kmag(p.ny, p.nx, p.dy, p.dx);
        build_active_modes(p, zeta);
        if (p.n_active <= 0) throw std::runtime_error("No active modes selected");
        build_interaction_plan(p, kmag);

        std::vector<cdouble> k1z(n), k1p(n), k2z(n), k2p(n), k3z(n), k3p(n), k4z(n), k4p(n);
        std::vector<cdouble> tmpz(n), tmpp(n);

        for (std::int64_t step = 0; step < p.n_steps; ++step) {
            rhs(p, zeta, phi, k1z, k1p);
            axpy(tmpz, zeta, 0.5 * h, k1z);
            axpy(tmpp, phi, 0.5 * h, k1p);

            rhs(p, tmpz, tmpp, k2z, k2p);
            axpy(tmpz, zeta, 0.5 * h, k2z);
            axpy(tmpp, phi, 0.5 * h, k2p);

            rhs(p, tmpz, tmpp, k3z, k3p);
            axpy(tmpz, zeta, h, k3z);
            axpy(tmpp, phi, h, k3p);

            rhs(p, tmpz, tmpp, k4z, k4p);

            for (std::size_t i = 0; i < n; ++i) {
                zeta[i] += (h / 6.0) * (k1z[i] + 2.0 * k2z[i] + 2.0 * k3z[i] + k4z[i]);
                phi[i] += (h / 6.0) * (k1p[i] + 2.0 * k2p[i] + 2.0 * k3p[i] + k4p[i]);
            }
        }

        write_complex_vector(dir / "zeta_final.bin", zeta);
        write_complex_vector(dir / "phi_final.bin", phi);

        const auto t1 = std::chrono::steady_clock::now();
        const double seconds = std::chrono::duration<double>(t1 - t0).count();
        std::ofstream summary(dir / "summary.txt");
        summary << "n_total " << p.n_total << "\n";
        summary << "n_active " << p.n_active << "\n";
        summary << "n_lambda_steps " << p.n_steps << "\n";
#ifdef CREAMER_HAVE_OPENMP
        summary << "openmp_threads " << omp_get_max_threads() << "\n";
#else
        summary << "openmp_threads 1\n";
#endif
        summary << "runtime_s " << seconds << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "creamer_flow error: " << e.what() << "\n";
        return 1;
    }
}
