#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#ifdef CREAMER_HAVE_OPENMP
#include <omp.h>
#endif

namespace fs = std::filesystem;
using cdouble = std::complex<double>;

struct Var { char kind; int mode; };
struct Term { Var v[4]; int n = 0; double coeff = 0.0; };
struct NormalVar { char kind; int mode; };
struct DerivativeBuckets {
    std::vector<Term> linear_linear;
    std::vector<Term> one_target_one_linear;
};
using DerivativeBucketList = std::vector<DerivativeBuckets>;
struct Stats {
    std::int64_t n_h3_terms = 0, n_w3_terms = 0, n_h4_terms = 0, n_k4_terms = 0;
    std::int64_t used_terms = 0, candidate_terms = 0, resonant_terms = 0, nonresonant_terms = 0;
};

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

void write_complex_vector(const fs::path& path, const std::vector<cdouble>& data) {
    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("Could not open output: " + path.string());
    for (const auto& z : data) {
        const double re = z.real(), im = z.imag();
        out.write(reinterpret_cast<const char*>(&re), sizeof(double));
        out.write(reinterpret_cast<const char*>(&im), sizeof(double));
    }
}

void write_double_vector(const fs::path& path, const std::vector<double>& data) {
    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("Could not open output: " + path.string());
    out.write(reinterpret_cast<const char*>(data.data()), static_cast<std::streamsize>(data.size() * sizeof(double)));
}

int mode_to_index(int m, int nx) {
    int r = m % nx;
    if (r < 0) r += nx;
    return r;
}

int fft_mode_number(int idx, int nx) {
    return (idx <= nx / 2) ? idx : idx - nx;
}

double kernel_d_finite(double theta1, double theta2, double theta3,
                       double k1sq, double k2sq, double k3sq, double g) {
    const double den = theta1 * (theta2 + theta3 - theta1)
                     + theta2 * (theta1 + theta3 - theta2)
                     + theta3 * (theta1 + theta2 - theta3);
    if (std::abs(den) < 1e-14 || theta1 == 0.0 || theta2 == 0.0 || theta3 == 0.0) return 0.0;
    const double num = k1sq * (theta1 * theta1 - (theta2 - theta3) * (theta2 - theta3))
                     + k2sq * (theta2 * theta2 - (theta1 - theta3) * (theta1 - theta3))
                     + k3sq * (theta3 * theta3 - (theta1 - theta2) * (theta1 - theta2))
                     - 2.0 * theta1 * theta2 * theta3 * (theta1 + theta2 + theta3);
    return num / (12.0 * g * den);
}

double kernel_b_finite(double theta1, double theta2, double theta3,
                       double k1sq, double k2sq, double k3sq) {
    const double den = theta1 * (theta2 + theta3 - theta1)
                     + theta2 * (theta1 + theta3 - theta2)
                     + theta3 * (theta1 + theta2 - theta3);
    if (std::abs(den) < 1e-14) return 0.0;
    const double num = theta3 * (theta3 * (theta1 + theta2) - theta1 * theta1 - theta2 * theta2)
                     + theta3 * (k1sq + k2sq - 2.0 * k3sq)
                     + (theta1 - theta2) * (k1sq - k2sq);
    return num / (2.0 * den);
}

struct Context {
    int nx = 0;
    double dx = 1.0, g = 9.81, h = 1.0;
    std::vector<int> mx;
    std::vector<double> k, theta;
    std::unordered_map<int, int> mode_pos;
    std::vector<int> mode_set, linear_modes, positive_modes, targets;
    std::unordered_set<int> linear_set, target_set, mode_set_lookup;
};

struct LinearState {
    std::vector<cdouble> z;
    std::vector<cdouble> p;
    std::vector<cdouble> a;
    std::vector<cdouble> b;
    std::vector<double> gamma;
};

double theta_of(const Context& c, int m) {
    if (std::abs(m) > c.nx / 2) return 0.0;
    return c.theta[static_cast<std::size_t>(mode_to_index(m, c.nx))];
}
double k_of(const Context& c, int m) {
    if (std::abs(m) > c.nx / 2) return 0.0;
    return c.k[static_cast<std::size_t>(mode_to_index(m, c.nx))];
}
bool has_mode(const Context& c, int m) { return c.mode_set_lookup.find(m) != c.mode_set_lookup.end(); }

void local_bd(const Context& c, int dest, int a, int b, double& D, double& Bz, double& Bphi) {
    const double td = theta_of(c, dest), ta = theta_of(c, a), tb = theta_of(c, b);
    const double kd2 = k_of(c, dest) * k_of(c, dest);
    const double ka2 = k_of(c, a) * k_of(c, a);
    const double kb2 = k_of(c, b) * k_of(c, b);
    D = kernel_d_finite(td, ta, tb, kd2, ka2, kb2, c.g);
    Bz = kernel_b_finite(ta, tb, td, ka2, kb2, kd2);
    Bphi = kernel_b_finite(td, ta, tb, kd2, ka2, kb2);
}

std::vector<int> unique_sorted(std::vector<int> v) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    return v;
}

Context build_context(const std::vector<cdouble>& zeta, int nx, double dx, double g, double h, int max_active) {
    Context c;
    c.nx = nx; c.dx = dx; c.g = g; c.h = h;
    const double L = nx * dx;
    c.mx.resize(nx); c.k.resize(nx); c.theta.resize(nx);
    for (int i = 0; i < nx; ++i) {
        c.mx[i] = fft_mode_number(i, nx);
        c.k[i] = (2.0 * M_PI / L) * static_cast<double>(c.mx[i]);
        const double ak = std::abs(c.k[i]);
        c.theta[i] = (ak == 0.0) ? 0.0 : ak * std::tanh(ak * h);
    }
    std::vector<std::pair<double, int>> energy;
    for (int i = 0; i < nx; ++i) {
        if (c.mx[i] > 0) energy.emplace_back(std::norm(zeta[i]), c.mx[i]);
    }
    std::stable_sort(energy.begin(), energy.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
    const int keep = std::min<int>(max_active, static_cast<int>(energy.size()));
    for (int i = 0; i < keep; ++i) c.positive_modes.push_back(energy[static_cast<std::size_t>(i)].second);
    std::sort(c.positive_modes.begin(), c.positive_modes.end());
    for (int m : c.positive_modes) { c.linear_modes.push_back(m); c.linear_modes.push_back(-m); }
    c.linear_modes = unique_sorted(c.linear_modes);
    c.linear_set.insert(c.linear_modes.begin(), c.linear_modes.end());

    std::vector<int> ms = c.linear_modes;
    for (int a : c.linear_modes) for (int b : c.linear_modes) ms.push_back(a + b);
    for (int a : c.linear_modes) for (int b : c.linear_modes) for (int d : c.linear_modes) ms.push_back(a + b + d);
    for (int m : ms) if (m != 0 && std::abs(m) <= nx / 2) c.mode_set.push_back(m);
    c.mode_set = unique_sorted(c.mode_set);
    c.mode_set_lookup.insert(c.mode_set.begin(), c.mode_set.end());
    for (std::size_t i = 0; i < c.mode_set.size(); ++i) c.mode_pos[c.mode_set[i]] = static_cast<int>(i);

    std::vector<int> targets;
    for (int a : c.positive_modes) for (int b : c.positive_modes) for (int d : c.positive_modes) {
        const int t = a + b + d;
        if (std::abs(t) <= nx / 2) targets.push_back(t);
    }
    c.targets = unique_sorted(targets);
    c.target_set.insert(c.targets.begin(), c.targets.end());
    return c;
}

void build_linear_state(const Context& c, const std::vector<cdouble>& zeta, LinearState& s) {
    const std::size_t n = c.mode_set.size();
    s.z.assign(n, 0.0);
    s.p.assign(n, 0.0);
    s.a.assign(n, 0.0);
    s.b.assign(n, 0.0);
    s.gamma.assign(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const int m = c.mode_set[i];
        const int idx = mode_to_index(m, c.nx);
        s.z[i] = zeta[static_cast<std::size_t>(idx)];
        const double th = theta_of(c, m);
        if (th == 0.0) {
            s.p[i] = 0.0;
            s.gamma[i] = 0.0;
        } else if (m > 0) {
            s.p[i] = cdouble(0.0, -1.0) * std::sqrt(c.g / th) * s.z[i];
            s.gamma[i] = std::sqrt(th / c.g);
        } else {
            s.p[i] = cdouble(0.0, 1.0) * std::sqrt(c.g / th) * s.z[i];
            s.gamma[i] = std::sqrt(th / c.g);
        }
        s.a[i] = s.z[i] + cdouble(0.0, 1.0) * s.gamma[i] * s.p[i];
        s.b[i] = s.z[i] - cdouble(0.0, 1.0) * s.gamma[i] * s.p[i];
    }
}

cdouble map_get(const Context& c, const std::vector<cdouble>& v, int m) {
    auto it = c.mode_pos.find(m);
    if (it == c.mode_pos.end()) return 0.0;
    return v[static_cast<std::size_t>(it->second)];
}

double gamma_get(const Context& c, const LinearState& s, int m) {
    auto it = c.mode_pos.find(m);
    if (it == c.mode_pos.end()) return 0.0;
    return s.gamma[static_cast<std::size_t>(it->second)];
}

void compute_h3_part(const Context& c, const std::vector<cdouble>& zlin, const std::vector<cdouble>& plin,
                     std::vector<cdouble>& z3_h3) {
    const std::size_t n = c.mode_set.size();
    std::vector<cdouble> z2(n, 0.0), p2(n, 0.0);
    for (std::size_t id = 0; id < n; ++id) {
        const int dest = c.mode_set[id];
        cdouble zs = 0.0, ps = 0.0;
        for (int a : c.mode_set) {
            const int b = dest - a;
            if (!has_mode(c, b)) continue;
            double D, Bz, Bphi; local_bd(c, dest, a, b, D, Bz, Bphi);
            const cdouble za = map_get(c, zlin, a), zb = map_get(c, zlin, b);
            const cdouble pa = map_get(c, plin, a), pb = map_get(c, plin, b);
            zs += 3.0 * D * pa * pb + Bz * za * zb;
            ps += -2.0 * Bphi * za * pb;
        }
        z2[id] = zs; p2[id] = ps;
    }
    z3_h3.assign(c.nx, 0.0);
    for (int dest : c.targets) {
        cdouble zs = 0.0;
        for (int a : c.mode_set) {
            const int b = dest - a;
            if (!has_mode(c, b)) continue;
            double D, Bz, Bphi; local_bd(c, dest, a, b, D, Bz, Bphi);
            zs += 3.0 * D * (map_get(c, p2, a) * map_get(c, plin, b) + map_get(c, plin, a) * map_get(c, p2, b))
                + Bz * (map_get(c, z2, a) * map_get(c, zlin, b) + map_get(c, zlin, a) * map_get(c, z2, b));
        }
        z3_h3[static_cast<std::size_t>(mode_to_index(dest, c.nx))] += 0.5 * zs;
    }
}

std::vector<Term> build_h3_terms(const Context& c) {
    std::vector<Term> out;
    out.reserve(c.mode_set.size() * c.mode_set.size());
    for (int a : c.mode_set) for (int b : c.mode_set) {
        const int cc = -(a + b);
        if (!has_mode(c, cc)) continue;
        const double coeff = 0.25 * (k_of(c, a) * k_of(c, a) + k_of(c, b) * k_of(c, b)
            - k_of(c, cc) * k_of(c, cc) - 2.0 * theta_of(c, a) * theta_of(c, b));
        if (coeff == 0.0) continue;
        Term t; t.n = 3; t.coeff = coeff; t.v[0] = {'z', cc}; t.v[1] = {'p', a}; t.v[2] = {'p', b};
        out.push_back(t);
    }
    return out;
}

std::vector<Term> build_w3_terms(const Context& c) {
    std::vector<Term> out;
    out.reserve(2 * c.mode_set.size() * c.mode_set.size());
    for (int a : c.mode_set) for (int b : c.mode_set) {
        const int cc = -(a + b);
        if (!has_mode(c, cc)) continue;
        double D, Bz, Bphi; local_bd(c, a, b, cc, D, Bz, Bphi);
        if (D != 0.0) { Term t; t.n = 3; t.coeff = D; t.v[0] = {'p', a}; t.v[1] = {'p', b}; t.v[2] = {'p', cc}; out.push_back(t); }
        if (Bphi != 0.0) { Term t; t.n = 3; t.coeff = Bphi; t.v[0] = {'z', a}; t.v[1] = {'z', b}; t.v[2] = {'p', cc}; out.push_back(t); }
    }
    return out;
}

std::vector<Term> derivative_terms(const std::vector<Term>& terms, Var target) {
    std::vector<Term> out;
    out.reserve(terms.size() / 16 + 8);
    for (const auto& t : terms) {
        int first = -1, count = 0;
        for (int i = 0; i < t.n; ++i) {
            if (t.v[i].kind == target.kind && t.v[i].mode == target.mode) {
                if (first < 0) first = i;
                ++count;
            }
        }
        if (count == 0) continue;
        Term d; d.n = t.n - 1; d.coeff = t.coeff * static_cast<double>(count);
        int q = 0;
        for (int i = 0; i < t.n; ++i) if (i != first) d.v[q++] = t.v[i];
        out.push_back(d);
    }
    return out;
}

bool can_contribute(const Context& c, const int modes[4]) {
    for (int target : c.targets) {
        for (int i = 0; i < 4; ++i) if (modes[i] == -target) {
            bool ok = true;
            for (int j = 0; j < 4; ++j) if (j != i && c.linear_set.find(modes[j]) == c.linear_set.end()) { ok = false; break; }
            if (ok) return true;
        }
    }
    return false;
}

bool is_target_mode(const Context& c, int m) {
    return (m < 0) && (c.target_set.find(-m) != c.target_set.end());
}

bool is_linear_mode(const Context& c, int m) {
    return c.linear_set.find(m) != c.linear_set.end();
}

std::vector<std::vector<Term>> precompute_derivatives(const std::vector<Term>& terms, const Context& c, char kind) {
    std::vector<std::vector<Term>> out(c.mode_set.size());
    for (const auto& t : terms) {
        for (std::size_t ik = 0; ik < c.mode_set.size(); ++ik) {
            const int mode = c.mode_set[ik];
            int first = -1, count = 0;
            for (int i = 0; i < t.n; ++i) {
                if (t.v[i].kind == kind && t.v[i].mode == mode) {
                    if (first < 0) first = i;
                    ++count;
                }
            }
            if (count == 0) continue;
            Term d; d.n = t.n - 1; d.coeff = t.coeff * static_cast<double>(count);
            int q = 0;
            for (int i = 0; i < t.n; ++i) if (i != first) d.v[q++] = t.v[i];
            out[ik].push_back(d);
        }
    }
    return out;
}

DerivativeBuckets bucket_derivative_terms(const Context& c, const std::vector<Term>& terms) {
    DerivativeBuckets buckets;
    buckets.linear_linear.reserve(terms.size());
    buckets.one_target_one_linear.reserve(terms.size());
    for (const auto& t : terms) {
        const int m0 = t.v[0].mode;
        const int m1 = t.v[1].mode;
        const bool lin0 = is_linear_mode(c, m0);
        const bool lin1 = is_linear_mode(c, m1);
        const bool tar0 = is_target_mode(c, m0);
        const bool tar1 = is_target_mode(c, m1);
        if (lin0 && lin1) {
            buckets.linear_linear.push_back(t);
        } else if ((tar0 && lin1) || (tar1 && lin0)) {
            buckets.one_target_one_linear.push_back(t);
        }
    }
    return buckets;
}

DerivativeBucketList precompute_buckets(const Context& c, const std::vector<std::vector<Term>>& derivs) {
    DerivativeBucketList out(derivs.size());
    for (std::size_t i = 0; i < derivs.size(); ++i) {
        out[i] = bucket_derivative_terms(c, derivs[i]);
    }
    return out;
}

bool linear_normal_term(const Context& c, const LinearState& s, const Term& term, int deriv_pos, char deriv_kind,
                        NormalVar nv[4], cdouble& coeff) {
    coeff = term.coeff;
    for (int i = 0; i < 4; ++i) {
        const Var v = term.v[i];
        char kind = deriv_kind;
        if (i != deriv_pos) {
            if (v.mode > 0) kind = 'a';
            else if (v.mode < 0) kind = 'b';
            else return false;
        }
        const double gamma = gamma_get(c, s, v.mode);
        cdouble piece = 0.0;
        if (v.kind == 'z') piece = 0.5;
        else if (gamma == 0.0) return false;
        else if (kind == 'a') piece = 1.0 / (cdouble(0.0, 2.0) * gamma);
        else piece = -1.0 / (cdouble(0.0, 2.0) * gamma);
        coeff *= piece;
        nv[i] = {kind, v.mode};
    }
    return true;
}

bool linear_normal_vars(const Context& c, const LinearState& s, const Var vars[4], double coeff_in,
                        int deriv_pos, char deriv_kind, NormalVar nv[4], cdouble& coeff) {
    coeff = coeff_in;
    for (int i = 0; i < 4; ++i) {
        const Var v = vars[i];
        char kind = deriv_kind;
        if (i != deriv_pos) {
            if (v.mode > 0) kind = 'a';
            else if (v.mode < 0) kind = 'b';
            else return false;
        }
        const double gamma = gamma_get(c, s, v.mode);
        cdouble piece = 0.0;
        if (v.kind == 'z') piece = 0.5;
        else if (gamma == 0.0) return false;
        else if (kind == 'a') piece = 1.0 / (cdouble(0.0, 2.0) * gamma);
        else piece = -1.0 / (cdouble(0.0, 2.0) * gamma);
        coeff *= piece;
        nv[i] = {kind, v.mode};
    }
    return true;
}

double normal_weight(const Context& c, const NormalVar nv[4]) {
    double w = 0.0;
    for (int i = 0; i < 4; ++i) {
        const double om = std::sqrt(c.g * theta_of(c, nv[i].mode));
        w += (nv[i].kind == 'a') ? om : -om;
    }
    return w;
}

cdouble eval_normal_var(const Context& c, const LinearState& s, const NormalVar& v) {
    auto it = c.mode_pos.find(v.mode);
    if (it == c.mode_pos.end()) return 0.0;
    const std::size_t idx = static_cast<std::size_t>(it->second);
    return (v.kind == 'a') ? s.a[idx] : s.b[idx];
}

void accumulate_w4_term(const Context& c, const LinearState& s, const Term& term,
                        std::vector<cdouble>& z3_h4, Stats& stats) {
    int modes[4];
    for (int i = 0; i < 4; ++i) modes[i] = term.v[i].mode;
    if (!can_contribute(c, modes)) return;
    for (int target : c.targets) {
        for (int deriv = 0; deriv < 4; ++deriv) if (modes[deriv] == -target) {
            bool ok = true;
            for (int j = 0; j < 4; ++j) if (j != deriv && c.linear_set.find(modes[j]) == c.linear_set.end()) { ok = false; break; }
            if (!ok) continue;
            for (char deriv_kind : {'a', 'b'}) {
                NormalVar nv[4]; cdouble normal_coeff;
                if (!linear_normal_term(c, s, term, deriv, deriv_kind, nv, normal_coeff)) continue;
                ++stats.candidate_terms;
                const double weight = normal_weight(c, nv);
                if (std::abs(weight) < 1e-10) { ++stats.resonant_terms; continue; }
                ++stats.nonresonant_terms; ++stats.used_terms;
                const cdouble wcoeff = -normal_coeff / (cdouble(0.0, 1.0) * weight);
                const double gamma = gamma_get(c, s, nv[deriv].mode);
                const cdouble dval = (nv[deriv].kind == 'a') ? cdouble(0.0, 1.0) * gamma : cdouble(0.0, -1.0) * gamma;
                cdouble val = wcoeff * dval;
                for (int j = 0; j < 4; ++j) if (j != deriv) val *= eval_normal_var(c, s, nv[j]);
                z3_h4[static_cast<std::size_t>(mode_to_index(target, c.nx))] += val;
            }
        }
    }
}

bool locate_single_target_term(const Context& c, const Term& term, int& deriv, int& target) {
    deriv = -1;
    target = 0;
    for (int i = 0; i < 4; ++i) {
        if (is_target_mode(c, term.v[i].mode)) {
            if (deriv >= 0) return false;
            deriv = i;
            target = -term.v[i].mode;
        } else if (!is_linear_mode(c, term.v[i].mode)) {
            return false;
        }
    }
    return deriv >= 0;
}

void accumulate_w4_term_known(const Context& c, const LinearState& s, const Term& term, int deriv, int target,
                              std::vector<cdouble>& z3_h4, Stats& stats) {
    for (char deriv_kind : {'a', 'b'}) {
        NormalVar nv[4]; cdouble normal_coeff;
        if (!linear_normal_term(c, s, term, deriv, deriv_kind, nv, normal_coeff)) continue;
        ++stats.candidate_terms;
        const double weight = normal_weight(c, nv);
        if (std::abs(weight) < 1e-10) { ++stats.resonant_terms; continue; }
        ++stats.nonresonant_terms; ++stats.used_terms;
        const cdouble wcoeff = -normal_coeff / (cdouble(0.0, 1.0) * weight);
        const double gamma = gamma_get(c, s, nv[deriv].mode);
        const cdouble dval = (nv[deriv].kind == 'a') ? cdouble(0.0, 1.0) * gamma : cdouble(0.0, -1.0) * gamma;
        cdouble val = wcoeff * dval;
        for (int j = 0; j < 4; ++j) if (j != deriv) val *= eval_normal_var(c, s, nv[j]);
        z3_h4[static_cast<std::size_t>(mode_to_index(target, c.nx))] += val;
    }
}

void accumulate_w4_term_fast(const Context& c, const LinearState& s, const Term& term,
                             std::vector<cdouble>& z3_h4, Stats& stats) {
    int deriv = -1;
    int target = 0;
    if (!locate_single_target_term(c, term, deriv, target)) return;
    accumulate_w4_term_known(c, s, term, deriv, target, z3_h4, stats);
}

void accumulate_pair_product_direct(const Context& c, const LinearState& s,
                                    const Term& ta, const Term& tb, double coeff_scale,
                                    bool target_in_a, std::vector<cdouble>& z3_h4, Stats& stats) {
    Var vars[4] = {ta.v[0], ta.v[1], tb.v[0], tb.v[1]};
    int deriv = -1;
    int target = 0;
    if (target_in_a) {
        if (is_target_mode(c, ta.v[0].mode)) { deriv = 0; target = -ta.v[0].mode; }
        else if (is_target_mode(c, ta.v[1].mode)) { deriv = 1; target = -ta.v[1].mode; }
        else return;
    } else {
        if (is_target_mode(c, tb.v[0].mode)) { deriv = 2; target = -tb.v[0].mode; }
        else if (is_target_mode(c, tb.v[1].mode)) { deriv = 3; target = -tb.v[1].mode; }
        else return;
    }

    const double coeff_raw = coeff_scale * ta.coeff * tb.coeff;
    for (char deriv_kind : {'a', 'b'}) {
        NormalVar nv[4]; cdouble normal_coeff;
        if (!linear_normal_vars(c, s, vars, coeff_raw, deriv, deriv_kind, nv, normal_coeff)) continue;
        ++stats.candidate_terms;
        const double weight = normal_weight(c, nv);
        if (std::abs(weight) < 1e-10) { ++stats.resonant_terms; continue; }
        ++stats.nonresonant_terms; ++stats.used_terms;
        const cdouble wcoeff = -normal_coeff / (cdouble(0.0, 1.0) * weight);
        const double gamma = gamma_get(c, s, nv[deriv].mode);
        const cdouble dval = (nv[deriv].kind == 'a') ? cdouble(0.0, 1.0) * gamma : cdouble(0.0, -1.0) * gamma;
        cdouble val = wcoeff * dval;
        for (int j = 0; j < 4; ++j) if (j != deriv) val *= eval_normal_var(c, s, nv[j]);
        z3_h4[static_cast<std::size_t>(mode_to_index(target, c.nx))] += val;
    }
}

void compute_h4_direct_part(const Context& c, const LinearState& s, std::vector<cdouble>& z3_h4, Stats& stats) {
    for (int target : c.targets) {
        // term: p[-target] z[b] z[c] p[d], with c = target - b - d
        for (int b : c.linear_modes) for (int d : c.linear_modes) {
            const int cc = target - b - d;
            if (!is_linear_mode(c, cc)) continue;
            Term t; t.n = 4;
            t.v[0] = {'p', -target}; t.v[1] = {'z', b}; t.v[2] = {'z', cc}; t.v[3] = {'p', d};
            t.coeff = 0.5 * (theta_of(c, target) * theta_of(c, target - b) * theta_of(c, d)
                - 0.5 * theta_of(c, target) * k_of(c, d) * k_of(c, d)
                - 0.5 * k_of(c, target) * k_of(c, target) * theta_of(c, d));
            if (t.coeff == 0.0) continue;
            ++stats.n_h4_terms; ++stats.n_k4_terms;
            accumulate_w4_term_known(c, s, t, 0, target, z3_h4, stats);
        }

        // term: p[a] z[-target] z[c] p[d], with c = -a + target - d
        for (int a : c.linear_modes) for (int d : c.linear_modes) {
            const int cc = -a + target - d;
            if (!is_linear_mode(c, cc)) continue;
            Term t; t.n = 4;
            t.v[0] = {'p', a}; t.v[1] = {'z', -target}; t.v[2] = {'z', cc}; t.v[3] = {'p', d};
            const int kk = -a;
            t.coeff = 0.5 * (theta_of(c, kk) * theta_of(c, kk + target) * theta_of(c, d)
                - 0.5 * theta_of(c, kk) * k_of(c, d) * k_of(c, d)
                - 0.5 * k_of(c, kk) * k_of(c, kk) * theta_of(c, d));
            if (t.coeff == 0.0) continue;
            ++stats.n_h4_terms; ++stats.n_k4_terms;
            accumulate_w4_term_known(c, s, t, 1, target, z3_h4, stats);
        }

        // term: p[a] z[b] z[-target] p[d], with d = target - a - b
        for (int a : c.linear_modes) for (int b : c.linear_modes) {
            const int d = target - a - b;
            if (!is_linear_mode(c, d)) continue;
            Term t; t.n = 4;
            t.v[0] = {'p', a}; t.v[1] = {'z', b}; t.v[2] = {'z', -target}; t.v[3] = {'p', d};
            const int kk = -a;
            t.coeff = 0.5 * (theta_of(c, kk) * theta_of(c, kk - b) * theta_of(c, d)
                - 0.5 * theta_of(c, kk) * k_of(c, d) * k_of(c, d)
                - 0.5 * k_of(c, kk) * k_of(c, kk) * theta_of(c, d));
            if (t.coeff == 0.0) continue;
            ++stats.n_h4_terms; ++stats.n_k4_terms;
            accumulate_w4_term_known(c, s, t, 2, target, z3_h4, stats);
        }

        // term: p[a] z[b] z[c] p[-target], with c = -a - b + target
        for (int a : c.linear_modes) for (int b : c.linear_modes) {
            const int cc = -a - b + target;
            if (!is_linear_mode(c, cc)) continue;
            Term t; t.n = 4;
            t.v[0] = {'p', a}; t.v[1] = {'z', b}; t.v[2] = {'z', cc}; t.v[3] = {'p', -target};
            const int kk = -a;
            t.coeff = 0.5 * (theta_of(c, kk) * theta_of(c, kk - b) * theta_of(c, -target)
                - 0.5 * theta_of(c, kk) * k_of(c, -target) * k_of(c, -target)
                - 0.5 * k_of(c, kk) * k_of(c, kk) * theta_of(c, -target));
            if (t.coeff == 0.0) continue;
            ++stats.n_h4_terms; ++stats.n_k4_terms;
            accumulate_w4_term_known(c, s, t, 3, target, z3_h4, stats);
        }
    }
}

void compute_h4_part(const Context& c, const LinearState& s,
                     std::vector<cdouble>& z3_h4, Stats& stats) {
    z3_h4.assign(c.nx, 0.0);
    compute_h4_direct_part(c, s, z3_h4, stats);

    const auto H3 = build_h3_terms(c);
    const auto W3 = build_w3_terms(c);
    stats.n_h3_terms = static_cast<std::int64_t>(H3.size());
    stats.n_w3_terms = static_cast<std::int64_t>(W3.size());
    const auto dH3z = precompute_derivatives(H3, c, 'z');
    const auto dH3p = precompute_derivatives(H3, c, 'p');
    const auto dW3z = precompute_derivatives(W3, c, 'z');
    const auto dW3p = precompute_derivatives(W3, c, 'p');
    const auto dW3z_buckets = precompute_buckets(c, dW3z);
    const auto dW3p_buckets = precompute_buckets(c, dW3p);
    std::vector<cdouble> bracket_z3(c.nx, 0.0);
    Stats bracket_stats;

#ifdef CREAMER_HAVE_OPENMP
#pragma omp parallel
    {
        std::vector<cdouble> local_z3(c.nx, 0.0);
        Stats local_stats;
#pragma omp for schedule(dynamic)
        for (int ik = 0; ik < static_cast<int>(c.mode_set.size()); ++ik) {
            const int k = c.mode_set[static_cast<std::size_t>(ik)];
#else
        for (int ik = 0; ik < static_cast<int>(c.mode_set.size()); ++ik) {
            const int k = c.mode_set[static_cast<std::size_t>(ik)];
            auto& local_z3 = bracket_z3;
            auto& local_stats = bracket_stats;
#endif
        const auto neg_it = c.mode_pos.find(-k);
        if (neg_it == c.mode_pos.end()) continue;
        const int ink = neg_it->second;
        for (int which = 0; which < 2; ++which) {
            const auto& A = (which == 0) ? dH3z[static_cast<std::size_t>(ik)] : dH3p[static_cast<std::size_t>(ik)];
            const auto& B = (which == 0) ? dW3p[static_cast<std::size_t>(ink)] : dW3z[static_cast<std::size_t>(ink)];
            const auto& buckets = (which == 0) ? dW3p_buckets[static_cast<std::size_t>(ink)] : dW3z_buckets[static_cast<std::size_t>(ink)];
            const double scale = (which == 0) ? -0.5 : 0.5;
            if (A.empty() || B.empty()) continue;
            for (const auto& ta : A) {
                const int ma[2] = {ta.v[0].mode, ta.v[1].mode};
                const bool a0_linear = is_linear_mode(c, ma[0]);
                const bool a1_linear = is_linear_mode(c, ma[1]);
                const bool a0_target = is_target_mode(c, ma[0]);
                const bool a1_target = is_target_mode(c, ma[1]);
                const std::vector<Term>* source_terms = nullptr;
                if (a0_linear && a1_linear) {
                    source_terms = &buckets.one_target_one_linear;
                } else if ((a0_target && a1_linear) || (a1_target && a0_linear)) {
                    source_terms = &buckets.linear_linear;
                } else {
                    continue;
                }
                const bool target_in_a = ((a0_target && a1_linear) || (a1_target && a0_linear));
                for (const auto& tb : *source_terms) {
                    ++local_stats.n_k4_terms;
                    accumulate_pair_product_direct(c, s, ta, tb, scale, target_in_a, local_z3, local_stats);
                }
            }
        }
    }
#ifdef CREAMER_HAVE_OPENMP
#pragma omp critical
        {
            for (int i = 0; i < c.nx; ++i) bracket_z3[static_cast<std::size_t>(i)] += local_z3[static_cast<std::size_t>(i)];
            bracket_stats.n_k4_terms += local_stats.n_k4_terms;
            bracket_stats.used_terms += local_stats.used_terms;
            bracket_stats.candidate_terms += local_stats.candidate_terms;
            bracket_stats.resonant_terms += local_stats.resonant_terms;
            bracket_stats.nonresonant_terms += local_stats.nonresonant_terms;
        }
    }
#endif
    for (int i = 0; i < c.nx; ++i) z3_h4[static_cast<std::size_t>(i)] += bracket_z3[static_cast<std::size_t>(i)];
    stats.n_k4_terms += bracket_stats.n_k4_terms;
    stats.used_terms += bracket_stats.used_terms;
    stats.candidate_terms += bracket_stats.candidate_terms;
    stats.resonant_terms += bracket_stats.resonant_terms;
    stats.nonresonant_terms += bracket_stats.nonresonant_terms;
}

int main(int argc, char** argv) {
    try {
        if (argc != 2) { std::cerr << "Usage: creamer_eta33_triad_1d <job_dir>\n"; return 2; }
        const auto t0 = std::chrono::steady_clock::now();
        const fs::path dir(argv[1]);
        const auto meta = read_vector<std::int64_t>(dir / "meta.bin", 2);
        const auto params = read_vector<double>(dir / "params.bin", 4);
        const int nx = static_cast<int>(meta[0]);
        const int max_active = static_cast<int>(meta[1]);
        const double dx = params[0], g = params[1], h = params[2];
        auto zeta = read_complex_vector(dir / "zeta0.bin", static_cast<std::size_t>(nx));

        Context c = build_context(zeta, nx, dx, g, h, max_active);
        LinearState s;
        std::vector<cdouble> eta33_h3, eta33_h4;
        build_linear_state(c, zeta, s);
        compute_h3_part(c, s.z, s.p, eta33_h3);
        Stats stats;
        compute_h4_part(c, s, eta33_h4, stats);

        for (int i = 1; i < nx / 2; ++i) {
            eta33_h3[static_cast<std::size_t>(nx - i)] = std::conj(eta33_h3[static_cast<std::size_t>(i)]);
            eta33_h4[static_cast<std::size_t>(nx - i)] = std::conj(eta33_h4[static_cast<std::size_t>(i)]);
        }
        eta33_h3[0] = eta33_h4[0] = 0.0;
        eta33_h3[static_cast<std::size_t>(nx / 2)] = cdouble(eta33_h3[static_cast<std::size_t>(nx / 2)].real(), 0.0);
        eta33_h4[static_cast<std::size_t>(nx / 2)] = cdouble(eta33_h4[static_cast<std::size_t>(nx / 2)].real(), 0.0);
        write_complex_vector(dir / "eta33_h3_hat.bin", eta33_h3);
        write_complex_vector(dir / "eta33_h4_delta_hat.bin", eta33_h4);
        write_double_vector(dir / "k.bin", c.k);
        write_double_vector(dir / "theta.bin", c.theta);

        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::ofstream summary(dir / "summary.txt");
        summary << "n_positive " << c.positive_modes.size() << "\n";
        summary << "n_modes " << c.mode_set.size() << "\n";
        summary << "n_targets " << c.targets.size() << "\n";
        summary << "n_h3_terms " << stats.n_h3_terms << "\n";
        summary << "n_w3_terms " << stats.n_w3_terms << "\n";
        summary << "n_h4_terms " << stats.n_h4_terms << "\n";
        summary << "n_k4_terms " << stats.n_k4_terms << "\n";
        summary << "used_terms " << stats.used_terms << "\n";
        summary << "candidate_terms " << stats.candidate_terms << "\n";
#ifdef CREAMER_HAVE_OPENMP
        summary << "openmp_threads " << omp_get_max_threads() << "\n";
#else
        summary << "openmp_threads 1\n";
#endif
        summary << "runtime_s " << seconds << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "creamer_eta33_triad_1d error: " << e.what() << "\n";
        return 1;
    }
}
