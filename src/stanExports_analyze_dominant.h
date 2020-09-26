// Generated by rstantools.  Do not edit by hand.

/*
    Hickory is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Hickory is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Hickory.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.21.0
#include <stan/model/model_header.hpp>
namespace model_analyze_dominant_namespace {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
static int current_statement_begin__;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_analyze_dominant");
    reader.add_event(80, 78, "end", "model_analyze_dominant");
    return reader;
}
#include <stan_meta_header.hpp>
class model_analyze_dominant
  : public stan::model::model_base_crtp<model_analyze_dominant> {
private:
        int N_loci;
        int N_pops;
        std::vector<std::vector<int> > n;
        std::vector<std::vector<int> > N;
        double mu_pi;
        double sd_pi;
        double mu_f;
        double sd_f;
        double mu_theta;
        double sd_theta;
public:
    model_analyze_dominant(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_analyze_dominant(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, random_seed__, pstream__);
    }
    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning
        current_statement_begin__ = -1;
        static const char* function__ = "model_analyze_dominant_namespace::model_analyze_dominant";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 6;
            context__.validate_dims("data initialization", "N_loci", "int", context__.to_vec());
            N_loci = int(0);
            vals_i__ = context__.vals_i("N_loci");
            pos__ = 0;
            N_loci = vals_i__[pos__++];
            check_greater_or_equal(function__, "N_loci", N_loci, 0);
            current_statement_begin__ = 7;
            context__.validate_dims("data initialization", "N_pops", "int", context__.to_vec());
            N_pops = int(0);
            vals_i__ = context__.vals_i("N_pops");
            pos__ = 0;
            N_pops = vals_i__[pos__++];
            check_greater_or_equal(function__, "N_pops", N_pops, 0);
            current_statement_begin__ = 8;
            validate_non_negative_index("n", "N_pops", N_pops);
            validate_non_negative_index("n", "N_loci", N_loci);
            context__.validate_dims("data initialization", "n", "int", context__.to_vec(N_pops,N_loci));
            n = std::vector<std::vector<int> >(N_pops, std::vector<int>(N_loci, int(0)));
            vals_i__ = context__.vals_i("n");
            pos__ = 0;
            size_t n_k_0_max__ = N_pops;
            size_t n_k_1_max__ = N_loci;
            for (size_t k_1__ = 0; k_1__ < n_k_1_max__; ++k_1__) {
                for (size_t k_0__ = 0; k_0__ < n_k_0_max__; ++k_0__) {
                    n[k_0__][k_1__] = vals_i__[pos__++];
                }
            }
            size_t n_i_0_max__ = N_pops;
            size_t n_i_1_max__ = N_loci;
            for (size_t i_0__ = 0; i_0__ < n_i_0_max__; ++i_0__) {
                for (size_t i_1__ = 0; i_1__ < n_i_1_max__; ++i_1__) {
                    check_greater_or_equal(function__, "n[i_0__][i_1__]", n[i_0__][i_1__], 0);
                }
            }
            current_statement_begin__ = 9;
            validate_non_negative_index("N", "N_pops", N_pops);
            validate_non_negative_index("N", "N_loci", N_loci);
            context__.validate_dims("data initialization", "N", "int", context__.to_vec(N_pops,N_loci));
            N = std::vector<std::vector<int> >(N_pops, std::vector<int>(N_loci, int(0)));
            vals_i__ = context__.vals_i("N");
            pos__ = 0;
            size_t N_k_0_max__ = N_pops;
            size_t N_k_1_max__ = N_loci;
            for (size_t k_1__ = 0; k_1__ < N_k_1_max__; ++k_1__) {
                for (size_t k_0__ = 0; k_0__ < N_k_0_max__; ++k_0__) {
                    N[k_0__][k_1__] = vals_i__[pos__++];
                }
            }
            size_t N_i_0_max__ = N_pops;
            size_t N_i_1_max__ = N_loci;
            for (size_t i_0__ = 0; i_0__ < N_i_0_max__; ++i_0__) {
                for (size_t i_1__ = 0; i_1__ < N_i_1_max__; ++i_1__) {
                    check_greater_or_equal(function__, "N[i_0__][i_1__]", N[i_0__][i_1__], 0);
                }
            }
            current_statement_begin__ = 12;
            context__.validate_dims("data initialization", "mu_pi", "double", context__.to_vec());
            mu_pi = double(0);
            vals_r__ = context__.vals_r("mu_pi");
            pos__ = 0;
            mu_pi = vals_r__[pos__++];
            current_statement_begin__ = 13;
            context__.validate_dims("data initialization", "sd_pi", "double", context__.to_vec());
            sd_pi = double(0);
            vals_r__ = context__.vals_r("sd_pi");
            pos__ = 0;
            sd_pi = vals_r__[pos__++];
            check_greater_or_equal(function__, "sd_pi", sd_pi, 0);
            current_statement_begin__ = 14;
            context__.validate_dims("data initialization", "mu_f", "double", context__.to_vec());
            mu_f = double(0);
            vals_r__ = context__.vals_r("mu_f");
            pos__ = 0;
            mu_f = vals_r__[pos__++];
            current_statement_begin__ = 15;
            context__.validate_dims("data initialization", "sd_f", "double", context__.to_vec());
            sd_f = double(0);
            vals_r__ = context__.vals_r("sd_f");
            pos__ = 0;
            sd_f = vals_r__[pos__++];
            check_greater_or_equal(function__, "sd_f", sd_f, 0);
            current_statement_begin__ = 16;
            context__.validate_dims("data initialization", "mu_theta", "double", context__.to_vec());
            mu_theta = double(0);
            vals_r__ = context__.vals_r("mu_theta");
            pos__ = 0;
            mu_theta = vals_r__[pos__++];
            current_statement_begin__ = 17;
            context__.validate_dims("data initialization", "sd_theta", "double", context__.to_vec());
            sd_theta = double(0);
            vals_r__ = context__.vals_r("sd_theta");
            pos__ = 0;
            sd_theta = vals_r__[pos__++];
            check_greater_or_equal(function__, "sd_theta", sd_theta, 0);
            // initialize transformed data variables
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 21;
            num_params_r__ += 1;
            current_statement_begin__ = 22;
            num_params_r__ += 1;
            current_statement_begin__ = 23;
            validate_non_negative_index("logit_pi", "N_loci", N_loci);
            num_params_r__ += N_loci;
            current_statement_begin__ = 26;
            validate_non_negative_index("p", "N_loci", N_loci);
            validate_non_negative_index("p", "N_pops", N_pops);
            num_params_r__ += (N_loci * N_pops);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_analyze_dominant() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        typedef double local_scalar_t__;
        stan::io::writer<double> writer__(params_r__, params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        current_statement_begin__ = 21;
        if (!(context__.contains_r("logit_f")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable logit_f missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("logit_f");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "logit_f", "double", context__.to_vec());
        double logit_f(0);
        logit_f = vals_r__[pos__++];
        try {
            writer__.scalar_unconstrain(logit_f);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable logit_f: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 22;
        if (!(context__.contains_r("logit_theta")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable logit_theta missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("logit_theta");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "logit_theta", "double", context__.to_vec());
        double logit_theta(0);
        logit_theta = vals_r__[pos__++];
        try {
            writer__.scalar_unconstrain(logit_theta);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable logit_theta: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 23;
        if (!(context__.contains_r("logit_pi")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable logit_pi missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("logit_pi");
        pos__ = 0U;
        validate_non_negative_index("logit_pi", "N_loci", N_loci);
        context__.validate_dims("parameter initialization", "logit_pi", "vector_d", context__.to_vec(N_loci));
        Eigen::Matrix<double, Eigen::Dynamic, 1> logit_pi(N_loci);
        size_t logit_pi_j_1_max__ = N_loci;
        for (size_t j_1__ = 0; j_1__ < logit_pi_j_1_max__; ++j_1__) {
            logit_pi(j_1__) = vals_r__[pos__++];
        }
        try {
            writer__.vector_unconstrain(logit_pi);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable logit_pi: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 26;
        if (!(context__.contains_r("p")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable p missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("p");
        pos__ = 0U;
        validate_non_negative_index("p", "N_loci", N_loci);
        validate_non_negative_index("p", "N_pops", N_pops);
        context__.validate_dims("parameter initialization", "p", "vector_d", context__.to_vec(N_pops,N_loci));
        std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1> > p(N_pops, Eigen::Matrix<double, Eigen::Dynamic, 1>(N_loci));
        size_t p_j_1_max__ = N_loci;
        size_t p_k_0_max__ = N_pops;
        for (size_t j_1__ = 0; j_1__ < p_j_1_max__; ++j_1__) {
            for (size_t k_0__ = 0; k_0__ < p_k_0_max__; ++k_0__) {
                p[k_0__](j_1__) = vals_r__[pos__++];
            }
        }
        size_t p_i_0_max__ = N_pops;
        for (size_t i_0__ = 0; i_0__ < p_i_0_max__; ++i_0__) {
            try {
                writer__.vector_lub_unconstrain(0, 1, p[i_0__]);
            } catch (const std::exception& e) {
                stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable p: ") + e.what()), current_statement_begin__, prog_reader__());
            }
        }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }
    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(std::vector<T__>& params_r__,
                 std::vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // dummy to suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
            // model parameters
            current_statement_begin__ = 21;
            local_scalar_t__ logit_f;
            (void) logit_f;  // dummy to suppress unused var warning
            if (jacobian__)
                logit_f = in__.scalar_constrain(lp__);
            else
                logit_f = in__.scalar_constrain();
            current_statement_begin__ = 22;
            local_scalar_t__ logit_theta;
            (void) logit_theta;  // dummy to suppress unused var warning
            if (jacobian__)
                logit_theta = in__.scalar_constrain(lp__);
            else
                logit_theta = in__.scalar_constrain();
            current_statement_begin__ = 23;
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> logit_pi;
            (void) logit_pi;  // dummy to suppress unused var warning
            if (jacobian__)
                logit_pi = in__.vector_constrain(N_loci, lp__);
            else
                logit_pi = in__.vector_constrain(N_loci);
            current_statement_begin__ = 26;
            std::vector<Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> > p;
            size_t p_d_0_max__ = N_pops;
            p.reserve(p_d_0_max__);
            for (size_t d_0__ = 0; d_0__ < p_d_0_max__; ++d_0__) {
                if (jacobian__)
                    p.push_back(in__.vector_lub_constrain(0, 1, N_loci, lp__));
                else
                    p.push_back(in__.vector_lub_constrain(0, 1, N_loci));
            }
            // transformed parameters
            current_statement_begin__ = 30;
            local_scalar_t__ f;
            (void) f;  // dummy to suppress unused var warning
            stan::math::initialize(f, DUMMY_VAR__);
            stan::math::fill(f, DUMMY_VAR__);
            current_statement_begin__ = 31;
            local_scalar_t__ theta;
            (void) theta;  // dummy to suppress unused var warning
            stan::math::initialize(theta, DUMMY_VAR__);
            stan::math::fill(theta, DUMMY_VAR__);
            current_statement_begin__ = 32;
            validate_non_negative_index("pi", "N_loci", N_loci);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> pi(N_loci);
            stan::math::initialize(pi, DUMMY_VAR__);
            stan::math::fill(pi, DUMMY_VAR__);
            current_statement_begin__ = 33;
            validate_non_negative_index("x", "N_loci", N_loci);
            validate_non_negative_index("x", "N_pops", N_pops);
            std::vector<std::vector<local_scalar_t__> > x(N_loci, std::vector<local_scalar_t__>(N_pops, local_scalar_t__(0)));
            stan::math::initialize(x, DUMMY_VAR__);
            stan::math::fill(x, DUMMY_VAR__);
            // transformed parameters block statements
            current_statement_begin__ = 35;
            stan::math::assign(f, inv_logit(logit_f));
            current_statement_begin__ = 36;
            stan::math::assign(theta, inv_logit(logit_theta));
            current_statement_begin__ = 37;
            stan::math::assign(pi, inv_logit(logit_pi));
            current_statement_begin__ = 39;
            for (int i = 1; i <= N_loci; ++i) {
                current_statement_begin__ = 40;
                for (int j = 1; j <= N_pops; ++j) {
                    current_statement_begin__ = 41;
                    stan::model::assign(x, 
                                stan::model::cons_list(stan::model::index_uni(i), stan::model::cons_list(stan::model::index_uni(j), stan::model::nil_index_list())), 
                                (((pow(get_base1(get_base1(p, j, "p", 1), i, "p", 2), 2) * (1.0 - f)) + (f * get_base1(get_base1(p, j, "p", 1), i, "p", 2))) + (((2.0 * get_base1(get_base1(p, j, "p", 1), i, "p", 2)) * (1.0 - get_base1(get_base1(p, j, "p", 1), i, "p", 2))) * (1 - f))), 
                                "assigning variable x");
                }
            }
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            current_statement_begin__ = 30;
            if (stan::math::is_uninitialized(f)) {
                std::stringstream msg__;
                msg__ << "Undefined transformed parameter: f";
                stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable f: ") + msg__.str()), current_statement_begin__, prog_reader__());
            }
            check_greater_or_equal(function__, "f", f, 0);
            check_less_or_equal(function__, "f", f, 1);
            current_statement_begin__ = 31;
            if (stan::math::is_uninitialized(theta)) {
                std::stringstream msg__;
                msg__ << "Undefined transformed parameter: theta";
                stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable theta: ") + msg__.str()), current_statement_begin__, prog_reader__());
            }
            check_greater_or_equal(function__, "theta", theta, 0);
            check_less_or_equal(function__, "theta", theta, 1);
            current_statement_begin__ = 32;
            size_t pi_j_1_max__ = N_loci;
            for (size_t j_1__ = 0; j_1__ < pi_j_1_max__; ++j_1__) {
                if (stan::math::is_uninitialized(pi(j_1__))) {
                    std::stringstream msg__;
                    msg__ << "Undefined transformed parameter: pi" << "(" << j_1__ << ")";
                    stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable pi: ") + msg__.str()), current_statement_begin__, prog_reader__());
                }
            }
            check_greater_or_equal(function__, "pi", pi, 0);
            check_less_or_equal(function__, "pi", pi, 1);
            current_statement_begin__ = 33;
            size_t x_k_0_max__ = N_loci;
            size_t x_k_1_max__ = N_pops;
            for (size_t k_0__ = 0; k_0__ < x_k_0_max__; ++k_0__) {
                for (size_t k_1__ = 0; k_1__ < x_k_1_max__; ++k_1__) {
                    if (stan::math::is_uninitialized(x[k_0__][k_1__])) {
                        std::stringstream msg__;
                        msg__ << "Undefined transformed parameter: x" << "[" << k_0__ << "]" << "[" << k_1__ << "]";
                        stan::lang::rethrow_located(std::runtime_error(std::string("Error initializing variable x: ") + msg__.str()), current_statement_begin__, prog_reader__());
                    }
                }
            }
            size_t x_i_0_max__ = N_loci;
            size_t x_i_1_max__ = N_pops;
            for (size_t i_0__ = 0; i_0__ < x_i_0_max__; ++i_0__) {
                for (size_t i_1__ = 0; i_1__ < x_i_1_max__; ++i_1__) {
                    check_greater_or_equal(function__, "x[i_0__][i_1__]", x[i_0__][i_1__], 0);
                    check_less_or_equal(function__, "x[i_0__][i_1__]", x[i_0__][i_1__], 1);
                }
            }
            // model body
            current_statement_begin__ = 50;
            for (int i = 1; i <= N_loci; ++i) {
                current_statement_begin__ = 51;
                for (int j = 1; j <= N_pops; ++j) {
                    current_statement_begin__ = 52;
                    lp_accum__.add(binomial_log<propto__>(get_base1(get_base1(n, j, "n", 1), i, "n", 2), get_base1(get_base1(N, j, "N", 1), i, "N", 2), get_base1(get_base1(x, i, "x", 1), j, "x", 2)));
                }
            }
            current_statement_begin__ = 58;
            lp_accum__.add(normal_log<propto__>(logit_pi, mu_pi, sd_pi));
            current_statement_begin__ = 59;
            lp_accum__.add(normal_log<propto__>(logit_f, mu_f, sd_f));
            current_statement_begin__ = 60;
            lp_accum__.add(normal_log<propto__>(logit_theta, mu_theta, sd_theta));
            current_statement_begin__ = 61;
            for (int i = 1; i <= N_loci; ++i) {
                current_statement_begin__ = 62;
                for (int j = 1; j <= N_pops; ++j) {
                    current_statement_begin__ = 63;
                    lp_accum__.add(beta_log<propto__>(get_base1(get_base1(p, j, "p", 1), i, "p", 2), (((1.0 - theta) / theta) * get_base1(pi, i, "pi", 1)), (((1.0 - theta) / theta) * (1.0 - get_base1(pi, i, "pi", 1)))));
                }
            }
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
        lp_accum__.add(lp__);
        return lp_accum__.sum();
    } // log_prob()
    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }
    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("logit_f");
        names__.push_back("logit_theta");
        names__.push_back("logit_pi");
        names__.push_back("p");
        names__.push_back("f");
        names__.push_back("theta");
        names__.push_back("pi");
        names__.push_back("x");
        names__.push_back("log_lik");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(N_loci);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(N_pops);
        dims__.push_back(N_loci);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(N_loci);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(N_loci);
        dims__.push_back(N_pops);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
    }
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;
        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
        static const char* function__ = "model_analyze_dominant_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        double logit_f = in__.scalar_constrain();
        vars__.push_back(logit_f);
        double logit_theta = in__.scalar_constrain();
        vars__.push_back(logit_theta);
        Eigen::Matrix<double, Eigen::Dynamic, 1> logit_pi = in__.vector_constrain(N_loci);
        size_t logit_pi_j_1_max__ = N_loci;
        for (size_t j_1__ = 0; j_1__ < logit_pi_j_1_max__; ++j_1__) {
            vars__.push_back(logit_pi(j_1__));
        }
        std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1> > p;
        size_t p_d_0_max__ = N_pops;
        p.reserve(p_d_0_max__);
        for (size_t d_0__ = 0; d_0__ < p_d_0_max__; ++d_0__) {
            p.push_back(in__.vector_lub_constrain(0, 1, N_loci));
        }
        size_t p_j_1_max__ = N_loci;
        size_t p_k_0_max__ = N_pops;
        for (size_t j_1__ = 0; j_1__ < p_j_1_max__; ++j_1__) {
            for (size_t k_0__ = 0; k_0__ < p_k_0_max__; ++k_0__) {
                vars__.push_back(p[k_0__](j_1__));
            }
        }
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            // declare and define transformed parameters
            current_statement_begin__ = 30;
            double f;
            (void) f;  // dummy to suppress unused var warning
            stan::math::initialize(f, DUMMY_VAR__);
            stan::math::fill(f, DUMMY_VAR__);
            current_statement_begin__ = 31;
            double theta;
            (void) theta;  // dummy to suppress unused var warning
            stan::math::initialize(theta, DUMMY_VAR__);
            stan::math::fill(theta, DUMMY_VAR__);
            current_statement_begin__ = 32;
            validate_non_negative_index("pi", "N_loci", N_loci);
            Eigen::Matrix<double, Eigen::Dynamic, 1> pi(N_loci);
            stan::math::initialize(pi, DUMMY_VAR__);
            stan::math::fill(pi, DUMMY_VAR__);
            current_statement_begin__ = 33;
            validate_non_negative_index("x", "N_loci", N_loci);
            validate_non_negative_index("x", "N_pops", N_pops);
            std::vector<std::vector<double> > x(N_loci, std::vector<double>(N_pops, double(0)));
            stan::math::initialize(x, DUMMY_VAR__);
            stan::math::fill(x, DUMMY_VAR__);
            // do transformed parameters statements
            current_statement_begin__ = 35;
            stan::math::assign(f, inv_logit(logit_f));
            current_statement_begin__ = 36;
            stan::math::assign(theta, inv_logit(logit_theta));
            current_statement_begin__ = 37;
            stan::math::assign(pi, inv_logit(logit_pi));
            current_statement_begin__ = 39;
            for (int i = 1; i <= N_loci; ++i) {
                current_statement_begin__ = 40;
                for (int j = 1; j <= N_pops; ++j) {
                    current_statement_begin__ = 41;
                    stan::model::assign(x, 
                                stan::model::cons_list(stan::model::index_uni(i), stan::model::cons_list(stan::model::index_uni(j), stan::model::nil_index_list())), 
                                (((pow(get_base1(get_base1(p, j, "p", 1), i, "p", 2), 2) * (1.0 - f)) + (f * get_base1(get_base1(p, j, "p", 1), i, "p", 2))) + (((2.0 * get_base1(get_base1(p, j, "p", 1), i, "p", 2)) * (1.0 - get_base1(get_base1(p, j, "p", 1), i, "p", 2))) * (1 - f))), 
                                "assigning variable x");
                }
            }
            if (!include_gqs__ && !include_tparams__) return;
            // validate transformed parameters
            const char* function__ = "validate transformed params";
            (void) function__;  // dummy to suppress unused var warning
            current_statement_begin__ = 30;
            check_greater_or_equal(function__, "f", f, 0);
            check_less_or_equal(function__, "f", f, 1);
            current_statement_begin__ = 31;
            check_greater_or_equal(function__, "theta", theta, 0);
            check_less_or_equal(function__, "theta", theta, 1);
            current_statement_begin__ = 32;
            check_greater_or_equal(function__, "pi", pi, 0);
            check_less_or_equal(function__, "pi", pi, 1);
            current_statement_begin__ = 33;
            size_t x_i_0_max__ = N_loci;
            size_t x_i_1_max__ = N_pops;
            for (size_t i_0__ = 0; i_0__ < x_i_0_max__; ++i_0__) {
                for (size_t i_1__ = 0; i_1__ < x_i_1_max__; ++i_1__) {
                    check_greater_or_equal(function__, "x[i_0__][i_1__]", x[i_0__][i_1__], 0);
                    check_less_or_equal(function__, "x[i_0__][i_1__]", x[i_0__][i_1__], 1);
                }
            }
            // write transformed parameters
            if (include_tparams__) {
                vars__.push_back(f);
                vars__.push_back(theta);
                size_t pi_j_1_max__ = N_loci;
                for (size_t j_1__ = 0; j_1__ < pi_j_1_max__; ++j_1__) {
                    vars__.push_back(pi(j_1__));
                }
                size_t x_k_0_max__ = N_loci;
                size_t x_k_1_max__ = N_pops;
                for (size_t k_1__ = 0; k_1__ < x_k_1_max__; ++k_1__) {
                    for (size_t k_0__ = 0; k_0__ < x_k_0_max__; ++k_0__) {
                        vars__.push_back(x[k_0__][k_1__]);
                    }
                }
            }
            if (!include_gqs__) return;
            // declare and define generated quantities
            current_statement_begin__ = 70;
            double log_lik;
            (void) log_lik;  // dummy to suppress unused var warning
            stan::math::initialize(log_lik, DUMMY_VAR__);
            stan::math::fill(log_lik, DUMMY_VAR__);
            // generated quantities statements
            current_statement_begin__ = 72;
            stan::math::assign(log_lik, 0.0);
            current_statement_begin__ = 73;
            for (int i = 1; i <= N_loci; ++i) {
                current_statement_begin__ = 74;
                for (int j = 1; j <= N_pops; ++j) {
                    current_statement_begin__ = 75;
                    stan::math::assign(log_lik, (log_lik + binomial_log(get_base1(get_base1(n, j, "n", 1), i, "n", 2), get_base1(get_base1(N, j, "N", 1), i, "N", 2), get_base1(get_base1(x, i, "x", 1), j, "x", 2))));
                }
            }
            // validate, write generated quantities
            current_statement_begin__ = 70;
            vars__.push_back(log_lik);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }
    std::string model_name() const {
        return "model_analyze_dominant";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "logit_f";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "logit_theta";
        param_names__.push_back(param_name_stream__.str());
        size_t logit_pi_j_1_max__ = N_loci;
        for (size_t j_1__ = 0; j_1__ < logit_pi_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "logit_pi" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t p_j_1_max__ = N_loci;
        size_t p_k_0_max__ = N_pops;
        for (size_t j_1__ = 0; j_1__ < p_j_1_max__; ++j_1__) {
            for (size_t k_0__ = 0; k_0__ < p_k_0_max__; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "p" << '.' << k_0__ + 1 << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "f";
            param_names__.push_back(param_name_stream__.str());
            param_name_stream__.str(std::string());
            param_name_stream__ << "theta";
            param_names__.push_back(param_name_stream__.str());
            size_t pi_j_1_max__ = N_loci;
            for (size_t j_1__ = 0; j_1__ < pi_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "pi" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
            size_t x_k_0_max__ = N_loci;
            size_t x_k_1_max__ = N_pops;
            for (size_t k_1__ = 0; k_1__ < x_k_1_max__; ++k_1__) {
                for (size_t k_0__ = 0; k_0__ < x_k_0_max__; ++k_0__) {
                    param_name_stream__.str(std::string());
                    param_name_stream__ << "x" << '.' << k_0__ + 1 << '.' << k_1__ + 1;
                    param_names__.push_back(param_name_stream__.str());
                }
            }
        }
        if (!include_gqs__) return;
        param_name_stream__.str(std::string());
        param_name_stream__ << "log_lik";
        param_names__.push_back(param_name_stream__.str());
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "logit_f";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "logit_theta";
        param_names__.push_back(param_name_stream__.str());
        size_t logit_pi_j_1_max__ = N_loci;
        for (size_t j_1__ = 0; j_1__ < logit_pi_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "logit_pi" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
        size_t p_j_1_max__ = N_loci;
        size_t p_k_0_max__ = N_pops;
        for (size_t j_1__ = 0; j_1__ < p_j_1_max__; ++j_1__) {
            for (size_t k_0__ = 0; k_0__ < p_k_0_max__; ++k_0__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "p" << '.' << k_0__ + 1 << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
        }
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "f";
            param_names__.push_back(param_name_stream__.str());
            param_name_stream__.str(std::string());
            param_name_stream__ << "theta";
            param_names__.push_back(param_name_stream__.str());
            size_t pi_j_1_max__ = N_loci;
            for (size_t j_1__ = 0; j_1__ < pi_j_1_max__; ++j_1__) {
                param_name_stream__.str(std::string());
                param_name_stream__ << "pi" << '.' << j_1__ + 1;
                param_names__.push_back(param_name_stream__.str());
            }
            size_t x_k_0_max__ = N_loci;
            size_t x_k_1_max__ = N_pops;
            for (size_t k_1__ = 0; k_1__ < x_k_1_max__; ++k_1__) {
                for (size_t k_0__ = 0; k_0__ < x_k_0_max__; ++k_0__) {
                    param_name_stream__.str(std::string());
                    param_name_stream__ << "x" << '.' << k_0__ + 1 << '.' << k_1__ + 1;
                    param_names__.push_back(param_name_stream__.str());
                }
            }
        }
        if (!include_gqs__) return;
        param_name_stream__.str(std::string());
        param_name_stream__ << "log_lik";
        param_names__.push_back(param_name_stream__.str());
    }
}; // model
}  // namespace
typedef model_analyze_dominant_namespace::model_analyze_dominant stan_model;
#ifndef USING_R
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
#endif
#endif
