'use strict';

import { create, all } from 'mathjs'
const math = create(all, {});

if (typeof(math) == 'undefined') {
    console.error("* math.js is not imported or pre-loaded, please check your environment before using.");
}

export const metajs = {

    tfxd6: function(x) {
        return math.divide(
            math.round(
                math.multiply(
                    x,
                    1000000
                )
            ),
            1000000
        );
    },

    backtransf: function(x, sm, value, n) {
        return x;
    },

    asin2p: function(x, n, value) {
        // the n is harmonic.mean 
        // if (sm == "IRFT") {
        //     if (is.metabind)
        //         harmonic.mean < - x$t.harmonic.mean.ma
        //     else
        //         harmonic.mean < - 1 / mean(1 / x$time)
        // } else {
        //     if (is.metabind)
        //         harmonic.mean < - x$n.harmonic.mean.ma
        //     else
        //         harmonic.mean < - 1 / mean(1 / x$n)
        // }
        // the `x$n` is the total of the event+total
            
        if (typeof(n) == 'undefined') {
            n = null;
        }
        if (typeof(value) == 'undefined') {
            value = 'mean';
        }

        var minimum = math.asin(0);
        var maximum = math.asin(1);

        if (n != null) {
            minimum = 0.5 * (Math.asin(Math.sqrt(0 / (n + 1))) + Math.asin(Math.sqrt((0 + 1) / (n + 1))))
            maximum = 0.5 * (Math.asin(Math.sqrt(n / (n + 1))) + Math.asin(Math.sqrt((n + 1) / (n + 1))))
        }

        if (n == null) {
            return math.dotPow(
                math.sin(x),
                2
            );
        } else {
            // 0.5 * (1 - sign(cos(2 * x[sel])) *
            //            sqrt(1 - ( sin(2 * x[sel]) +
            //                         ( sin(2 * x[sel]) - 1 / sin(2 * x[sel])) / n[sel]
            //                             
            //                     )^2
            //            )
            //       )
            return math.dotMultiply(
                0.5,
                math.subtract(
                    1,
                    math.dotMultiply(
                        math.sign(math.cos(math.dotMultiply(2, x))),
                        math.sqrt(
                            math.subtract(
                                1,
                                math.dotPow(
                                    math.add(
                                        math.sin(math.dotMultiply(2, x)),
                                        math.dotDivide(
                                            math.subtract(
                                                math.sin(math.dotMultiply(2, x)),
                                                math.dotDivide(
                                                    1,
                                                    math.sin(math.dotMultiply(2, x))
                                                )
                                            ),
                                            n
                                        )
                                    ),
                                    2
                                )
                            )
                        )
                    )
                )
            )
        }
    },

    expit: function(v) {
        return math.add(
            math.dotMultiply(
                0.5,
                math.tanh(
                    math.dotMultiply(
                        0.5,
                        v
                    )
                )
            ),
            0.5
        );
    },

    harmonic_mean: function(n) {
        return math.dotDivide(
            1,
            math.mean(
                math.dotDivide(1, n)
            )
        );
    },
    
    inv_calc: function(X, W) {
        return math.inv(
            math.multiply(
                math.multiply(
                    math.transpose(X), 
                    W
                ), 
                X
            )
        );
    },

    /**
     * Calculate I-Squared
     * 
     * For more information, check:
     * https://rdrr.io/cran/meta/src/R/meta-het.R
     */
    isquared: function(Q, df, level) {
        if (typeof(level) == 'undefined') {
            level = 0.95;
        }
        var tres = this.calcH(Q, df, level);

        var func_t = function(x) {
            return (x**2 - 1) / (x**2);
        }

        return {
            TE: func_t(tres.TE),
            lower: func_t(tres.lower),
            upper: func_t(tres.upper),
        }
    },

    /**
     * Calculate H
     * 
     * For more information, check:
     * https://rdrr.io/cran/meta/src/R/meta-het.R
     */
    calcH: function(Q, df, level) {
        if (typeof(level) == 'undefined') {
            level = 0.95;
        }
        //
        // Calculate H
        // Higgins & Thompson (2002), Statistics in Medicine, 21, 1539-58
        //
        var k = df + 1;
        var H = NaN;
        var selogH = NaN;

        if (isNaN(k)) {
            
        } else {
            if (k > 1) {
                H = math.sqrt(Q / (k - 1));
            } else {
                H = NaN;   
            }

            if (Q > k) {
                if (k >= 2) {
                    selogH = 0.5 * (math.log(Q) - math.log(k-1)) / (math.sqrt(2*Q) - math.sqrt(2*k - 3));
                } else {
                    selogH = NaN;
                }
            } else {
                if (k > 2) {
                    selogH = math.sqrt(1/(2*(k-2)) * (1 - 1/(3*(k-2)**2)));
                } else {
                    selogH = NaN;
                }
            }
            
        }

        var tres = this.ci95(
            math.log(math.max(
                H, 1
            )),
            selogH
        );

        return {
            TE: math.exp(tres.TE),
            lower: math.max(math.exp(tres.lower), 1),
            upper: math.max(math.exp(tres.upper), 1)
        }
    },
    
    /**
     * Calculation of confidence intervals on ci level of 0.95
     * 
     * For more information,
     * https://rdrr.io/cran/meta/src/R/ci.R
     */
    ci95: function(TE, seTE) {
        var level = 0.95;

        var lower = math.subtract(
            TE,
            math.multiply(
                1.96,
                seTE
            )
        );

        var upper = math.add(
            TE,
            math.multiply(
                1.96,
                seTE
            )
        );

        var statistic = math.dotDivide(TE, seTE);

        var pval = 2 * pnorm.compute(
            math.abs(
                statistic
            ), 
            false
        );

        return {
            TE: TE,
            seTE: seTE,
            pval: pval,
            level: level,
            lower: this.tfxd6(lower),
            upper: this.tfxd6(upper),
            statistic: statistic
        };
    },
    
    /**
     * Network Meta-analysis based Frequentist Method
     * 
     * Due to the complixty of input format, `netmeta` accepts list of dict object
     * [
     *     {sm: 1.2, lower: 1.1, upper: 1.3, 
     *      treat1: 'Nivo', treat2: 'Suni', 
     *      study: 'TREX', year: 2020}, 
     * ]
     * 
     * 
     * @param {Object} rs dataset
     * @param {Object} params configuration
     */
    netmeta: function(rs, params) {
        
        ///////////////////////////////////////////////////
        // (1) Check data and set arguments
        ///////////////////////////////////////////////////
        // Yes, you can use default settings
        if (typeof(params)=='undefined') {
            params = {};
        }
        
        if (!params.hasOwnProperty('rs_format')) {
            params['rs_format'] = 'dict';
        }
        
        if (!params.hasOwnProperty('reference_group')) {
            params['reference_group'] = '';
        }
        
        if (!params.hasOwnProperty('sm')) {
            params['sm'] = 'HR';
        }
        
        if (!params.hasOwnProperty('small_values')) {
            params['small_values'] = 'good';
        }

        if (!params.hasOwnProperty('method_tau')) {
            params['method.tau'] = 'DL';
        }

        // check if any treat1 == treat2

        // check subgraph
        
        
        ///////////////////////////////////////////////////
        // (2) Read data
        ///////////////////////////////////////////////////

        ///////////////////////////////////////////////////
        // (3) Store dataset
        ///////////////////////////////////////////////////
        var ds = {
            studlab: [],

            // pre data
            TE: [],
            seTE: [],
            treat1: [],
            treat2: [],
            
            // raw data
            event1: [],
            n1: [],
            event2: [],
            n2: [],

            sd1: [],
            sd2: [],

            time1: [],
            time2: []
        };

        ///////////////////////////////////////////////////
        // (4) Additional checks
        ///////////////////////////////////////////////////


        var ret = {};

        return ret;
    },

    netconnection: function(rs) {

    },

    /**
     * Meta-analysis of single proportions
     * 
     * The input `rs` is a list that contains the records.
     * 
     * [
     *     [event, n],
     *     ...
     * ]
     * 
     * Each records contains two values:
     *     - event: the number of events
     *     - n: the number of observations
     * 
     * For more information, check the R code
     * https://rdrr.io/cran/meta/src/R/metaprop.R
     */
    metaprop: function(rs, params) {
        
        ///////////////////////////////////////////////////
        // (1) Check and set arguments
        ///////////////////////////////////////////////////
        // Yes, you can use default settings
        if (typeof(params)=='undefined') {
            params = {};
        }
        
        if (!params.hasOwnProperty('rs_format')) {
            params['rs_format'] = 'list';
        }

        if (!params.hasOwnProperty('input_format')) {
            params['input_format'] = 'PRIM_CAT_RAW';
        }

        if (!params.hasOwnProperty('sm')) {
            params['sm'] = 'PFT';
        }

        var incr = 0.5;
        if (!params.hasOwnProperty('incr')) {
            params['incr'] = incr;
        } else {
            // check float
            incr = params['incr'];
        }

        if (!params.hasOwnProperty('incr_event')) {
            params['incr_event'] = params['incr'];
        }

        if (!params.hasOwnProperty('method')) {
            params['method'] = 'Inverse';
        }

        if (!params.hasOwnProperty('method_ci')) {
            params['method_ci'] = 'CP';
        }

        if (!params.hasOwnProperty('method_tau')) {
            params['method.tau'] = 'DL';
        }

        ///////////////////////////////////////////////////
        // (2) Read data
        ///////////////////////////////////////////////////
        var _rs = [];
        if (params['rs_format'] == 'dict') {
            // convert the dict format to list
            var _rs = [];
            for (let i = 0; i < rs.length; i++) {
                const r = rs[i];
                var _e = NaN;
                var _n = NaN;

                // the event and total can be defined in several ways
                if (r.hasOwnProperty('Et')) { _e = r.Et; }
                if (r.hasOwnProperty('event')) { _e = r.event; }
                if (r.hasOwnProperty('Nt')) { _n = r.Nt; }
                if (r.hasOwnProperty('total')) { _n = r.total; }
                if (r.hasOwnProperty('n')) { _n = r.n; }
                
                _rs.push([
                    _e, _n
                ]);
            }
        } else {
            _rs = rs;
        }


        ///////////////////////////////////////////////////
        // (3) Check length of variables
        ///////////////////////////////////////////////////


        ///////////////////////////////////////////////////
        // (4) Subset, exclude studies, and subgroup
        ///////////////////////////////////////////////////


        ///////////////////////////////////////////////////
        // (5) Store dataset
        ///////////////////////////////////////////////////
        var ds = {
            e: [],
            n: [],
            incr_e: []
        };
        // a mapping from ds index to rs index
        var d2r = {};
        for (let i = 0; i < _rs.length; i++) {
            const r = _rs[i];
            // check the r
            if (r[0] > r[1]) {
                // it's not possible that event > n
                continue;

            } else if (r[0] == 0) {
                // zero event??? increase both
                ds.incr_e.push( incr );

            } else {
                // for most case
                ds.incr_e.push( 0 );
            }

            // then save both
            ds.e.push( r[0] );
            ds.n.push( r[1] );

            // add this mapping
            d2r[ds.e.length - 1] = i;
        }

        // double check the length of records
        if (ds.e.length == 0) {
            // what??? nothing to do with 0 studies
            return this.__mk_metaprop_nan(rs, params);

        } else if (ds.e.length == 1) {
            // what??? only one study?
            // ok ...
        } else {
            // ok, more than 1
        }

        // add harmonic_mean to ds
        ds.harmonic_mean = this.harmonic_mean(ds.n);

        ///////////////////////////////////////////////////
        // (6) Subset analysis
        ///////////////////////////////////////////////////


        ///////////////////////////////////////////////////
        // (7) Calculate results for each study
        ///////////////////////////////////////////////////
        var TE = null;
        var seTE = null;

        if (params.sm == 'PLOGIT') {
            // TE <- log((event + incr.event) / (n - event + incr.event))
            // seTE <- sqrt(1 / (event + incr.event) +
            //             1 / ((n - event + incr.event)))
            TE = math.log(
                math.dotDivide(
                    math.add(ds.e, ds.incr_e),
                    math.add(
                        math.subtract(ds.n, ds.e),
                        ds.incr_e
                    )
                )
            );

            seTE = math.sqrt(
                math.add(
                    math.dotDivide(
                        1, 
                        math.add(ds.e, ds.incr_e)
                    ),
                    math.dotDivide(
                        1,
                        math.add(
                            math.subtract(ds.n, ds.e),
                            ds.incr_e
                        )
                    )
                )
            );

        } else if (params.sm == 'PFT') {
            // TE <- 0.5 * (asin(sqrt(event / (n + 1))) + asin(sqrt((event + 1) / (n + 1))))
            // seTE <- sqrt(1 / (4 * n + 2))
            // transf.null.effect <- asin(sqrt(null.effect))
            TE = math.dotMultiply(
                0.5,
                math.add(
                    math.asin(math.sqrt(
                        math.dotDivide(
                            ds.e, 
                            math.add(ds.n, 1)
                    ))),
                    math.asin(math.sqrt(
                        math.dotDivide(
                            math.add(ds.e, 1), 
                            math.add(ds.n, 1)
                    )))
                )
            );
            seTE = math.sqrt(
                math.dotDivide(
                    1,
                    math.add(
                        2,
                        math.dotMultiply(
                            4,
                            ds.n
                        )
                    )
                )
            );

        } else if (params.sm == 'PAS') {
            // TE <- asin(sqrt(event / n))
            // seTE <- sqrt(1 / (4 * n))
            // transf.null.effect <- asin(sqrt(null.effect))

            TE = math.asin(math.sqrt(
                math.dotDivide(
                    ds.e, ds.n
                )
            ));

            seTE = math.sqrt(
                math.dotDivide(
                    1,
                    math.dotMultiply(4, ds.n)
                )
            );

        } else if (params.sm == 'PLN') {
            // TE <- log((event + incr.event) / (n + incr.event))
            //  Hartung, Knapp (2001), p. 3880, formula (18):
            // seTE <- ifelse(event == n,
            //             sqrt(1 / event                - 1 / (n + incr.event)),
            //             sqrt(1 / (event + incr.event) - 1 / (n + incr.event))
            //             )
            TE = math.log(
                math.dotDivide(
                    math.add(ds.e, ds.incr_e),
                    math.add(ds.n, ds.incr_e)
                )
            );

            if (math.deepEqual(ds.e, ds.n)) {
                // event == total??
                seTE = math.sqrt(
                    math.subtract(
                        math.dotDivide(1, ds.e),
                        math.dotDivide(
                            1,
                            math.add(ds.n, ds.incr_e)
                        )
                    )
                );
            } else {
                seTE = math.sqrt(
                    math.subtract(
                        math.dotDivide(
                            1, 
                            math.add(ds.e, ds.incr_e)
                        ),
                        math.dotDivide(
                            1,
                            math.add(ds.n, ds.incr_e)
                        )
                    )
                )
            }
        }

        var SM = [];
        var SM_lower = [];
        var SM_upper = [];

        // var SM_RS = ds.e.tolist().map((e,i)=>binom.test(e, ds.n.get(i)));
        for (let i = 0; i < ds.e.length; i++) {
            var _e = ds.e[i];
            var _incr_e = ds.incr_e[i];
            var _n = ds.n[i];
            var r = binom.test(_e + _incr_e, _n);

            SM.push(r.estimate);
            SM_lower.push(r.lower);
            SM_upper.push(r.upper);
        }

        // var _TE = Array(rs.length).fill(null);
        // var _seTE = Array(rs.length).fill(null);

        ///////////////////////////////////////////////////
        // (8) Do meta-analysis
        ///////////////////////////////////////////////////
        var fixed = null;
        var random = null;

        ///////////////////////////////
        // Fixed effect model
        ///////////////////////////////
        // var w_fixed = ONE.divide(seTE.pow(2));
        // var wp_fixed = w_fixed.divide(w_fixed.sum());
        var w_fixed = math.dotDivide(
            1, 
            math.dotPow(seTE, 2)
        );
        var wp_fixed = math.dotDivide(
            w_fixed,
            math.sum(w_fixed)
        );

        // var TE_fixed = TE.multiply(wp_fixed).sum();
        // var seTE_fixed = math.sqrt( 1 / w_fixed.sum());
        var TE_fixed = math.sum(
            math.dotMultiply(
                TE,
                wp_fixed
            )
        );
        var seTE_fixed = math.sqrt(
            math.divide(
                1,
                math.sum(w_fixed)
            )
        );
        var TE_fixed_lower = math.add(
            TE_fixed,
            math.dotMultiply(
                -1.96,
                seTE_fixed
            )
        );
        var TE_fixed_upper = math.add(
            TE_fixed,
            math.dotMultiply(
                1.96,
                seTE_fixed
            )
        );
        
        // var SM_fixed = this.expit(TE_fixed);
        // var SM_fixed_lower = this.expit();
        // var SM_fixed_upper = this.expit(TE_fixed + 1.96 * seTE_fixed);
        var SM_fixed = null;
        var SM_fixed_lower = null;
        var SM_fixed_upper = null;

        if (params.sm == 'PLOGIT') {
            SM_fixed = this.expit(TE_fixed);
            SM_fixed_lower = this.expit(TE_fixed_lower);
            SM_fixed_upper = this.expit(TE_fixed_upper);

        } else if (params.sm == 'PFT') {
            SM_fixed = this.asin2p(TE_fixed, ds.harmonic_mean);
            SM_fixed_lower = this.asin2p(TE_fixed_lower, ds.harmonic_mean);
            SM_fixed_upper = this.asin2p(TE_fixed_upper, ds.harmonic_mean);
        }

        fixed = {
            TE: TE_fixed,
            seTE: seTE_fixed,
            TE_lower: TE_fixed_lower,
            TE_upper: TE_fixed_upper,
            w: w_fixed,
            wp: wp_fixed,

            SM: SM_fixed,
            SM_lower: SM_fixed_lower,
            SM_upper: SM_fixed_upper
        }

        ///////////////////////////////
        // Heteroginity
        ///////////////////////////////
        var het = null;
        // het = this.__calc_heterogeneity_by_DL_with_TE_tau(TE, seTE, fixed.TE);
        het = this.__calc_heterogeneity_by_DL_rma_uni(TE, seTE);

        ///////////////////////////////
        // Random effect model
        ///////////////////////////////

        // for TE.length == 1, just use same as fixed
        if (TE.length == 1) {
            // for only one study, random is same 
            random = JSON.parse(JSON.stringify(fixed));
        } else {
            if (params.sm == 'PLOGIT') {
                random = this.__calc_random(TE, seTE, het);

                // backtransf
                random.SM = this.expit(random.TE);
                random.SM_lower = this.expit(random.TE_lower);
                random.SM_upper = this.expit(random.TE_upper);

            } else if (params.sm == 'PFT') {
                random = this.__calc_random(TE, seTE, het);

                // backtransf
                random.SM = this.asin2p(random.TE, ds.harmonic_mean);
                random.SM_lower = this.asin2p(random.TE_lower, ds.harmonic_mean);
                random.SM_upper = this.asin2p(random.TE_upper, ds.harmonic_mean);
            }
        }

        ///////////////////////////////////////////////////
        // (9) Finalize return object
        ///////////////////////////////////////////////////
        var ret = {
            ds: ds,
            heterogeneity: het,

            TE: TE,
            seTE: seTE,
            SM: SM,
            SM_lower: SM_lower,
            SM_upper: SM_upper,

            // the MA result
            fixed: fixed,
            random: random,

            // the settings
            params: params
        }
        ret = this.expand_ret(ret, rs, d2r);

        return ret;
    },

    /**
     * Meta-analysis of binary outcome data
     * 
     * The input `rs` is a list that contains the records.
     * 
     * [
     *     [Et, Nt, Ec, Nc],
     *     ...
     * ]
     * 
     * Each records contains two values:
     *     - Et: the number of events in the treatment group
     *     - Nt: the number of observations in the treatment group
     *     - Ec: the number of events in the control group
     *     - Nc: the number of observations in the control group
     * 
     * For more information, check the R code
     * https://rdrr.io/cran/meta/src/R/metabin.R
     */
    metabin: function(rs, params) {
        ///////////////////////////////////////////////////
        // (1) Check and set arguments
        ///////////////////////////////////////////////////
        // Yes, you can use default settings
        if (typeof(params)=='undefined') {
            params = {};
        }
        
        if (!params.hasOwnProperty('rs_format')) {
            params['rs_format'] = 'list';
        }

        if (!params.hasOwnProperty('input_format')) {
            params['input_format'] = 'PRIM_CAT_RAW';
        }

        if (!params.hasOwnProperty('sm')) {
            params['sm'] = 'OR';
        }

        var incr = 0.5;
        if (!params.hasOwnProperty('incr')) {
            params['incr'] = incr;
        } else {
            // check float
            incr = params['incr'];
        }

        if (!params.hasOwnProperty('incr_event')) {
            params['incr_event'] = params['incr'];
        }

        if (!params.hasOwnProperty('method')) {
            params['method'] = 'MH';
        }

        if (!params.hasOwnProperty('method_bias')) {
            params['method_bias'] = 'Harbord';
        }

        if (!params.hasOwnProperty('method_ci')) {
            params['method_ci'] = 'CP';
        }

        if (!params.hasOwnProperty('method_tau')) {
            params['method_tau'] = 'DL';
        }

        if (!params.hasOwnProperty('method_tau_ci')) {
            params['method_tau_ci'] = 'J';
        }

        if (!params.hasOwnProperty('prediction')) {
            params['prediction'] = false;
        }

        if (!params.hasOwnProperty('model_glmm')) {
            params['model_glmm'] = 'UM.FS';
        }

        ///////////////////////////////////////////////////
        // (2) Read data
        ///////////////////////////////////////////////////
        var _rs = [];
        if (params['rs_format'] == 'dict') {
            // convert the dict format to list
            var _rs = [];
            for (let i = 0; i < rs.length; i++) {
                const r = rs[i];
                _rs.push([
                    r.Et, r.Nt, r.Ec, r.Nc
                ]);
            }
        } else {
            _rs = rs;
        }

        ///////////////////////////////////////////////////
        // (3) Check length of variables
        ///////////////////////////////////////////////////


        ///////////////////////////////////////////////////
        // (4) Subset, exclude studies, and subgroup
        ///////////////////////////////////////////////////


        ///////////////////////////////////////////////////
        // (5) Store dataset
        ///////////////////////////////////////////////////
        var ds = {
            Et: [],
            Nt: [],
            Ec: [],
            Nc: [],
            incr_e: [],
            incr_c: []
        };
        // a mapping from ds index to rs index
        var d2r = {};
        for (let i = 0; i < _rs.length; i++) {
            const r = _rs[i];
            // check the r
            if (r[0] > r[1]) {
                // it's not possible that event > n
                continue;

            } else if (r[0] == 0 && r[2] == 0) {
                // both 0???
                continue;

            } else if (params.sm == 'OR' && r[0] == r[1] && r[2] == r[3]) {
                // both all events???
                continue;

            } else if (r[0] == 0 || r[2] == 0) {
                // zero event??? increase both
                ds.incr_e.push( incr );
                ds.incr_c.push( incr );

            } else {
                // for most case, no incr
                ds.incr_e.push( 0 );
                ds.incr_c.push( 0 );
            }

            // now save the e and n
            ds.Et.push( r[0] );
            ds.Nt.push( r[1] );
            ds.Ec.push( r[2] );
            ds.Nc.push( r[3] );

            // add this mapping
            d2r[ds.Et.length - 1] = i;
        }
        // double check the length of records
        if (ds.Et.length == 0) {
            // what??? no qualified records???
            return this.__mk_metaprop_nan(rs, params);

        } else if (ds.Et.length == 1) {
            // what??? only one study?
            // well, still workable

        } else {
            // ok, more than 1
            // let's do math!
        }

        // for the calculation
        ds.n11 = ds.Et;
        ds.n21 = ds.Ec;
        ds.n1_ = ds.Nt;
        ds.n2_ = ds.Nc;

        // other for calculation
        ds.n__ = math.add(ds.Nt, ds.Nc);
        ds.n12 = math.subtract(ds.Nt, ds.Et);
        ds.n22 = math.subtract(ds.Nc, ds.Ec);
        ds.n_1 = math.add(ds.n11, ds.n21);
        ds.n_2 = math.add(ds.n12, ds.n22);

        ///////////////////////////////////////////////////
        // (6) Subset analysis
        ///////////////////////////////////////////////////


        ///////////////////////////////////////////////////
        // (7) Calculate results for each study
        ///////////////////////////////////////////////////
        var TE = null;
        var seTE = null;
        var SM = null;
        var SM_lower = null;
        var SM_upper = null;

        if (params.sm == 'OR') {
            if (['MH', 'Inverse', 'GLMM', 'SSW'].includes(params.method)) {
                // Cooper & Hedges (1994), p. 251-2
                TE = math.log(
                    math.dotDivide(
                        math.dotMultiply(
                            math.add(ds.n11, ds.incr_e),
                            math.add(ds.n22, ds.incr_c)
                        ),
                        math.dotMultiply(
                            math.add(ds.n12, ds.incr_e),
                            math.add(ds.n21, ds.incr_c)
                        ),
                    )
                );

                seTE = math.sqrt(
                    math.add(
                        math.dotDivide(1, math.add(ds.n11, ds.incr_e)),
                        math.dotDivide(1, math.add(ds.n12, ds.incr_e)),
                        math.dotDivide(1, math.add(ds.n21, ds.incr_c)),
                        math.dotDivide(1, math.add(ds.n22, ds.incr_c))
                    )
                );

                // backtransf
                SM = math.exp(TE);
                SM_lower = math.exp(
                    math.subtract(
                        TE,
                        math.multiply(
                            1.96,
                            seTE
                        )
                    )
                );
                SM_upper = math.exp(
                    math.add(
                        TE,
                        math.multiply(
                            1.96,
                            seTE
                        )
                    )
                );
            }

        } else if (params.sm == 'RR') {
            // TE = nj.log(
            //     (ds.n11.divide(ds.n1_)).divide(
            //         ds.n21.divide(ds.n2_)
            //     )
            // );
            TE = math.log(
                math.dotDivide(
                    math.dotDivide(
                        math.add(ds.n11, ds.incr_e),
                        math.add(ds.n1_, ds.incr_e)
                    ),
                    math.dotDivide(
                        math.add(ds.n21, ds.incr_c),
                        math.add(ds.n2_, ds.incr_c)
                    ),
                )
            );
            // Hartung & Knapp (2001), Stat Med, equation (18)
            // seTE = (
            //     ONE.divide(ds.n11)
            //     .subtract(ONE.divide(ds.n1_))
            //     .add(ONE.divide(ds.n21))
            //     .subtract(ONE.divide(ds.n2_))
            // ).pow(0.5);
            seTE = math.sqrt(
                math.add(
                    math.dotDivide(
                        1,
                        math.add(ds.n11, ds.incr_e)
                    ),
                    math.dotDivide(
                        -1,
                        math.add(ds.n1_, ds.incr_e)
                    ),
                    math.dotDivide(
                        1,
                        math.add(ds.n21, ds.incr_c)
                    ),
                    math.dotDivide(
                        -1,
                        math.add(ds.n2_, ds.incr_c)
                    )
                )
            );

            // backtransf 
            // SM = nj.exp(TE);
            // SM_lower = nj.exp(TE.subtract(seTE.multiply(1.96)));
            // SM_upper = nj.exp(TE.add(seTE.multiply(1.96)));
            SM = math.exp(TE);
            SM_lower = math.exp(
                math.subtract(
                    TE,
                    math.dotMultiply(
                        1.96,
                        seTE
                    )
                )
            );
            SM_upper = math.exp(
                math.add(
                    TE,
                    math.dotMultiply(
                        1.96,
                        seTE
                    )
                )
            );
        }

        var _TE = Array(rs.length).fill(null);
        var _seTE = Array(rs.length).fill(null);

        ///////////////////////////////////////////////////
        // (8) Do meta-analysis
        ///////////////////////////////////////////////////
        var fixed = null;
        var random = null;

        ///////////////////////////////
        // Fixed effect model
        ///////////////////////////////
        if (['MH'].includes(params.method)) {
            if (params.sm == 'OR') {
                fixed = this.__calc_fixed_OR_by_MH(ds);

            } else if (params.sm == 'RR') {
                fixed = this.__calc_fixed_RR_by_MH(ds);
            }
        }
        
        ///////////////////////////////
        // Heteroginity
        ///////////////////////////////
        var het = null;
        // heterogeneity = this.__calc_heterogeneity_by_DL(TE, seTE);
        het = this.__calc_heterogeneity_by_DL_with_TE_tau(
            TE, 
            seTE, 
            fixed.TE
        );

        ///////////////////////////////
        // Random effect model
        ///////////////////////////////

        if (params.method == 'SSW') {
            // metagen.R
            // w.random <- n.e * n.c / (n.e + n.c)
            // w.random[exclude] <- 0
            // TE.random <- weighted.mean(TE, w.random, na.rm = TRUE)
            // seTE.random <- sqrt(sum(w.random^2 * (seTE^2 + m$tau^2), na.rm = TRUE) /
            //                     sum(w.random, na.rm = TRUE)^2)
            // ##
            // w.random[is.na(w.random)] <- 0

        } else {
            // metagen.R
            // w.random <- 1/(seTE^2 + sum(tau2.calc, na.rm = TRUE))
            // w.random[is.na(w.random) | is.na(TE) | exclude] <- 0
            // TE.random <- weighted.mean(TE, w.random, na.rm = TRUE)
            // seTE.random <- sqrt(1/sum(w.random, na.rm = TRUE))
            if (TE.length == 1) {
                // for only one study, random is same 
                random = JSON.parse(JSON.stringify(fixed));

            } else {
                // most cases
                random = this.__calc_random(TE, seTE, het);

                // backtransf
                random.SM = math.exp(random.TE);
                random.SM_lower = math.exp(random.TE_lower);
                random.SM_upper = math.exp(random.TE_upper);
            }
        }

        ///////////////////////////////////////////////////
        // (9) Finalize return object
        ///////////////////////////////////////////////////
        var ret = {
            ds: ds,
            heterogeneity: het,

            // each study
            TE: TE,
            seTE: seTE,
            SM: SM,
            SM_lower: SM_lower,
            SM_upper: SM_upper,

            // the MA result
            fixed: fixed,
            random: random,

            // the settings
            params: params
        }

        // expand the ret to it's original size of given
        ret = this.expand_ret(ret, rs, d2r);

        return ret;
    },

    __mk_heterogeneity_nan: function() {
        return {
            I2: NaN, tau2: NaN, pval_Q: NaN
        }
    },

    __mk_metaprop_nan: function(rs, params) {
        var NaNs = Array(rs.length).fill(NaN);
        var _params = JSON.parse(JSON.stringify(params));

        return {
            ds: {},
            heterogeneity: this.__mk_heterogeneity_nan(),

            // each study
            TE: NaNs,
            seTE: NaNs,
            SM: NaNs,
            SM_lower: NaNs,
            SM_upper: NaNs,

            // the MA result
            fixed: {
                TE: NaN,
                seTE: NaN,
                TE_lower: NaN,
                TE_upper: NaN,
                w: NaNs,
                wp: NaNs,

                SM: NaN,
                SM_lower: NaN,
                SM_upper: NaN
            },
            random: {
                TE: NaN,
                seTE: NaN,
                TE_lower: NaN,
                TE_upper: NaN,
                w: NaNs,
                wp: NaNs,

                SM: NaN,
                SM_lower: NaN,
                SM_upper: NaN
            },

            // the settings
            params: _params
        }
    },

    __calc_fixed_OR_by_MH: function(ds) {
        var A = math.dotDivide(
            math.dotMultiply(
                math.add(ds.n11, ds.incr_e),
                math.add(ds.n22, ds.incr_c)
            ),
            math.add(
                ds.n__,
                math.dotMultiply(2, ds.incr_e),
                math.dotMultiply(2, ds.incr_c)
            )
        );

        var B = math.dotDivide(
            math.add(
                ds.n11,
                ds.incr_e,
                ds.n22,
                ds.incr_c
            ),
            math.add(
                ds.n__,
                math.dotMultiply(2, ds.incr_e),
                math.dotMultiply(2, ds.incr_c)
            )
        );

        var C = math.dotDivide(
            math.dotMultiply(
                math.add(ds.n12, ds.incr_e),
                math.add(ds.n21, ds.incr_c)
            ),
            math.add(
                ds.n__,
                math.dotMultiply(2, ds.incr_e),
                math.dotMultiply(2, ds.incr_c)
            )
        );

        var D = math.dotDivide(
            math.add(
                ds.n12,
                ds.incr_e,
                ds.n21,
                ds.incr_c
            ),
            math.add(
                ds.n__,
                math.dotMultiply(2, ds.incr_e),
                math.dotMultiply(2, ds.incr_c)
            )
        );

        var w_fixed = C;
        var wp_fixed = math.dotDivide(
            w_fixed,
            math.sum(w_fixed)
        );

        // TODO remove NaN from sum
        var TE_fixed = math.log(
            math.divide(
                math.sum(A),
                math.sum(C)
            )
        );

        var seTE_fixed = math.sqrt(
            (1 / (2 * (math.sum(A))**2)) * (
                math.sum(math.dotMultiply(A, B)) 
                +
                math.exp(TE_fixed) * (
                    math.sum(math.dotMultiply(B, C)) +
                    math.sum(math.dotMultiply(A, D))
                ) 
                +
                math.exp(TE_fixed)**2 * math.sum(math.dotMultiply(C, D))
            )
        );

        var TE_fixed_lower = TE_fixed - 1.96 * seTE_fixed;
        var TE_fixed_upper = TE_fixed + 1.96 * seTE_fixed

        var SM_fixed = math.exp(TE_fixed);
        var SM_fixed_lower = math.exp(TE_fixed_lower);
        var SM_fixed_upper = math.exp(TE_fixed_upper);

        var fixed = {
            TE: TE_fixed,
            seTE: seTE_fixed,
            TE_lower: TE_fixed_lower,
            TE_upper: TE_fixed_upper,
            w: w_fixed,
            wp: wp_fixed,

            SM: SM_fixed,
            SM_lower: SM_fixed_lower,
            SM_upper: SM_fixed_upper
        };

        return fixed;
    },


    __calc_fixed_RR_by_MH: function(ds) {
        var D = math.dotDivide(
            math.subtract(
                math.dotMultiply(
                    math.add(
                        ds.n1_,
                        math.dotMultiply(2, ds.incr_e)
                    ),
                    math.dotMultiply(
                        math.add(
                            ds.n2_, 
                            math.dotMultiply(2, ds.incr_c)
                        ),
                        math.add(
                            ds.n_1,
                            ds.incr_e, 
                            ds.incr_c
                        )
                    )
                ),
                math.dotMultiply(
                    math.add(ds.n11, ds.incr_e),
                    math.dotMultiply(
                        math.add(ds.n21, ds.incr_c),
                        math.add(
                            ds.n__,
                            math.dotMultiply(2, ds.incr_e),
                            math.dotMultiply(2, ds.incr_c)
                        )
                    )
                )
            ),
            math.dotPow(
                math.add(
                    ds.n__,
                    math.dotMultiply(2, ds.incr_e),
                    math.dotMultiply(2, ds.incr_c)
                ), 
                2
            )
        );

        var R = math.dotDivide(
            math.dotMultiply(
                math.add(ds.n11, ds.incr_e),
                math.add(ds.n2_, math.dotMultiply(2, ds.incr_c))
            ),
            math.add(
                ds.n__,
                math.dotMultiply(2, ds.incr_e),
                math.dotMultiply(2, ds.incr_c)
            )
        );

        var S = math.dotDivide(
            math.dotMultiply(
                math.add(ds.n21, ds.incr_c),
                math.add(ds.n1_, math.dotMultiply(2, ds.incr_e))
            ),
            math.add(
                ds.n__,
                math.dotMultiply(2, ds.incr_e),
                math.dotMultiply(2, ds.incr_c)
            )
        );

        var w_fixed = S;
        var wp_fixed = math.dotDivide(
            w_fixed,
            math.sum(w_fixed)
        );

        var TE_fixed = math.log(
            math.divide(
                math.sum(R),
                math.sum(S)
            )
        );

        var seTE_fixed = math.sqrt(
            math.divide(
                math.sum(D),
                math.multiply(
                    math.sum(R),
                    math.sum(S)
                )
            )
        );

        var TE_fixed_lower = TE_fixed - 1.96 * seTE_fixed;
        var TE_fixed_upper = TE_fixed + 1.96 * seTE_fixed

        var SM_fixed = math.exp(TE_fixed);
        var SM_fixed_lower = math.exp(TE_fixed_lower);
        var SM_fixed_upper = math.exp(TE_fixed_upper);

        var fixed = {
            TE: TE_fixed,
            seTE: seTE_fixed,
            TE_lower: TE_fixed_lower,
            TE_upper: TE_fixed_upper,
            w: w_fixed,
            wp: wp_fixed,

            SM: SM_fixed,
            SM_lower: SM_fixed_lower,
            SM_upper: SM_fixed_upper
        };

        return fixed;
    },


    __calc_random: function(TE, seTE, het) {
        var w_random = math.dotDivide(
            1,
            math.add(
                math.sum(het.tau2),
                math.dotPow(seTE, 2)
            )
        );
        w_random = this.fillna(w_random);

        var wp_random = math.dotDivide(
            w_random,
            math.sum(w_random)
        );

        var TE_random = math.sum(
            math.dotMultiply(TE, wp_random)
        );

        var seTE_random = math.sqrt(
            math.dotDivide(
                1,
                math.sum(w_random)
            )
        );

        var TE_random_lower = math.add(
            TE_random,
            math.dotMultiply(
                -1.96,
                seTE_random
            )
        );
        var TE_random_upper = math.add(
            TE_random,
            math.dotMultiply(
                1.96,
                seTE_random
            )
        );

        var random = {
            TE: TE_random,
            seTE: seTE_random,
            TE_lower: TE_random_lower,
            TE_upper: TE_random_upper,
            w: w_random,
            wp: wp_random,

            SM: null,
            SM_lower: null,
            SM_upper: null
        };

        return random;
    },

    /**
     * Heterogeneity estimation for standard model (rma.uni)
     * by the DerSimonian-Laird (DL) estimator
     * 
     * For more information:
     * https://rdrr.io/cran/meta/src/R/hetcalc.R
     * and 
     * https://rdrr.io/cran/metafor/src/R/rma.uni.r
     */
    __calc_heterogeneity_by_DL_rma_uni: function(TE, seTE) {
        // wi <- 1/vi
        // W <- diag(wi, nrow = k, ncol = k)
        // stXWX <- .invcalc(X = X, W = W, k = k)
        // P <- W - W %*% X %*% stXWX %*% crossprod(X, W)
        // RSS <- crossprod(Y, P) %*% Y
        // trP <- .tr(P)
        // tau2 <- ifelse(tau2.fix, tau2.val, (RSS - (k - p))/trP)

        if (TE.length == 1) {
            // if this is only one study, no heter
            return this.__mk_heterogeneity_nan();
        }

        var k = TE.length;
        var p = 1;

        var vi = math.dotPow(seTE, 2);
        var wi = math.dotDivide(1, vi);
        var W = math.diag(wi);
        
        var X = math.ones(TE.length, 1);
        var Y = math.transpose([TE]);
        var stXWX = this.inv_calc(X, W);

        var P = math.subtract(
            W, 
            math.multiply(
                math.multiply(
                    math.multiply(W, X),
                    stXWX
                ),
                math.multiply(
                    math.transpose(X),
                    W
                )
            )
        );
        var RSS = math.multiply(
            math.multiply(
                math.transpose(Y), 
                P
            ), 
            Y
        );
        RSS = RSS.toArray()[0][0];

        // the 
        var QE = math.max(0, RSS);
        var df_Q = k - p;

        var pval_Q = pchisq.compute(QE, df_Q, false);
        var I2 = this.isquared(QE, df_Q)
        
        var trP = math.sum(math.diag(P));
        var tau2 = (RSS - (k - p)) / trP;

        // make sure that tau2 is >= tau2_min
        // avoid neg val
        var tau2_min = 0;
        tau2 = math.max(tau2_min, tau2);

        // se.tau2 <- sqrt(1/trP^2 * (2*(k-p) + 4*max(tau2,0)*trP + 2*max(tau2,0)^2*sum(P*P)))
        var se_tau2 = math.sqrt(
            1 / trP**2 * (
                2 * (k - p) + 
                4 * math.max(tau2, 0) * trP +
                2 * math.max(tau2, 0) ** 2 * math.sum(math.multiply(P, P))
            )
        );

        // not sure about this:
        // I2 <- 100 * mean(tau2)/(vt + mean(tau2))

        var heterogeneity = {
            I2: I2.TE,
            tau2: tau2,
            pval_Q: pval_Q
        };

        return heterogeneity;
    },

    /**
     * Heterogeneity estimation
     * by the DerSimonian-Laird (DL) estimator
     * 
     * For more information:
     * https://rdrr.io/cran/meta/src/R/hetcalc.R
     */
    __calc_heterogeneity_by_DL_with_TE_tau: function(TE, seTE, TE_tau) {

        // Mantel-Haenszel estimator to calculate Q and tau (like RevMan 5)
        // w.fixed <- 1 / seTE^2
        // w.fixed[is.na(w.fixed)] <- 0
        // 
        // Q <- sum(w.fixed * (TE - TE.tau)^2, na.rm = TRUE)
        // df.Q <- sum(!is.na(seTE)) - 1
        // pval.Q <- pvalQ(Q, df.Q)
        // 
        // if (df.Q == 0)
        // tau2 <- NA
        // else if (round(Q, digits = 18) <= df.Q)
        // tau2 <- 0
        // else
        // tau2 <- (Q - df.Q) / Ccalc(w.fixed)

        if (TE.length == 1) {
            // if this is only one study, no heter
            return this.__mk_heterogeneity_nan();
        }

        function Ccalc(x) {
            return math.sum(x) - 
            math.sum(math.dotPow(x, 2)) / math.sum(x)
        }

        var w_fixed = math.dotDivide(
            1,
            math.dotPow(
                seTE,
                2
            )
        );

        // TE_tau is TE_fixed

        // fill na to 0
        w_fixed = this.fillna(w_fixed);

        // get Q
        var Q = math.sum(
            math.dotMultiply(
                w_fixed,
                math.dotPow(
                    math.subtract(TE, TE_tau),
                    2
                )
            )
        );

        var df_Q = math.subtract(
            math.sum(
                math.not(
                    math.isNaN(seTE)
                )
            ),
            1
        );

        var pval_Q = pchisq.compute(Q, df_Q, false);

        var tau2 = NaN;
        if (df_Q == 0) {
            tau2 = NaN;
        } else if (Q <= df_Q) {
            tau2 = 0;
        } else {
            tau2 = (Q - df_Q) / Ccalc(w_fixed)
        }

        // var tau = math.sqrt(tau2);

        // var H = this.calcH(Q, df_Q);

        var I2 = this.isquared(Q, df_Q);

        return {
            I2: I2.TE,
            tau2: tau2,
            pval_Q: pval_Q
        }
    },

    fillna: function(x) {
        if (typeof(x) == 'number') {
            return math.isNaN(x)?0:x;
        }
        return x.map(v=>math.isNaN(v)?0:v);
    },

    expand_list: function(vals, n, d2r) {
        var vs = new Array(n).fill(NaN);
        for (let i = 0; i < vals.length; i++) {
            vs[d2r[i]] = vs[i];
        }
        return vs;
    },

    expand_ret: function(ret, rs, d2r) {
        if (rs.length == Object.values(d2r).length) {
            // which means the records are mapped to ret.ds 1:1
            return ret;
        }
        // expand basic items
        var attrs = ['TE', 'seTE', 'SM', 'SM_lower', 'SM_upper'];
        for (let i = 0; i < attrs.length; i++) {
            ret[attrs[i]] = this.expand_list(ret[attrs[i]], rs.length, d2r);
        }
        // expand fixed and random weights
        ret.fixed.w = this.expand_list(ret.fixed.w, rs.length, d2r);
        ret.fixed.wp = this.expand_list(ret.fixed.wp, rs.length, d2r);
        ret.random.w = this.expand_list(ret.random.w, rs.length, d2r);
        ret.random.wp = this.expand_list(ret.random.wp, rs.length, d2r);

        return ret;
    }
}

export const pnorm = {

    normal_cdf: function(X) {
        var T = math.dotDivide(
            1,
            math.add(
                1,
                math.multiply(
                    0.2316419,
                    math.abs(X)
                )
            )
        );
        
        // var D=.3989423*Math.exp(-X*X/2);
        
        var D = math.multiply(
            0.3989423,
            math.exp(
                math.multiply(
                    -0.5,
                    math.dotPow(X, 2)
                )
            )
        );

        // var p=D*T*(.3193815+T*(-.3565638+T*(1.781478+T*(-1.821256+T*1.330274))));

        var p = math.dotMultiply(
            D,
            math.dotMultiply(
                T,
                math.add(
                    .3193815,
                    math.dotMultiply(
                        T,
                        math.add(
                            -.3565638,
                            math.dotMultiply(
                                T,
                                math.add(
                                    1.781478,
                                    math.dotMultiply(
                                        T,
                                        math.add(
                                            -1.821256,
                                            math.multiply(
                                                1.330274,
                                                T
                                            )
                                        )
                                    )
                                )
                            )
                        )
                    )
                )
            )
        );
        // if (X>0) {
        //     p = 1 - p;
        // }

        // we just know the X is always above 0
        p = math.subtract(
            1,
            p
        );
        return p
    },

    compute: function(q, lower_tail) {
        // mean is 0
        // sd is 1
        // otherwise
        // p = this.normal_cdf((q-mean)/sd)
        var p = this.normal_cdf(q);

        if (lower_tail) {

        } else {
            p = math.subtract(
                1, p
            )
        }

        p = math.round(100000 * p) / 100000;

        return p;
    }
};

export const pchisq = {
    // https://www.math.ucla.edu/~tom/distributions/chisq.html
    log_gamma: function(Z) {
        
        var S=1+76.18009173/Z-86.50532033/(Z+1)+24.01409822/(Z+2)-1.231739516/(Z+3)+.00120858003/(Z+4)-.00000536382/(Z+5);
        var LG= (Z-.5)*math.log(Z+4.5)-(Z+4.5)+math.log(S*2.50662827465);
    
        return LG
    },

    Gcf: function(X,A) {        
        // for X>A+1
        var A0=0;
        var B0=1;
        var A1=1;
        var B1=X;
        var AOLD=0;
        var N=0;
        while (math.abs((A1-AOLD)/A1)>.00001) {
            AOLD=A1;
            N=N+1;
            A0=A1+(N-A)*A0;
            B0=B1+(N-A)*B0;
            A1=X*A0+N*A1;
            B1=X*B0+N*B1;
            A0=A0/B1;
            B0=B0/B1;
            A1=A1/B1;
            B1=1;
        }
        var Prob=math.exp(A*math.log(X)-X- this.log_gamma(A))*A1;
    
        return 1-Prob;
    },

    Gser: function(X,A) {        
        // for X<A+1.
        var T9=1/A;
        var G=T9;
        var I=1;
        while (T9>G*.00001) {
            T9=T9*X/(A+I);
            G=G+T9;
            I=I+1;
        }
        G = G*math.exp(A*math.log(X)-X- this.log_gamma(A));

        return G;
    },

    gamma_cdf: function(x,a) {
        var GI;
        if (x<=0) {
            GI = 0;
        } else if (x<a+1) {
            GI = this.Gser(x,a);
        } else {
            GI = this.Gcf(x,a);
        }
        return GI;
    },

    compute: function(q, df, lower_tail) {
        if (typeof(lower_tail) == 'undefined') {
            lower_tail = true;
        }
        if (df<=0) {
            return NaN;
        }
        
        var p = this.gamma_cdf(q/2, df/2);
        p = math.round(p*100000) / 100000;

        if (lower_tail) {

        } else {
            p = 1 - p;
            p = math.round(p*100000) / 100000;
        }
        
        return p;
    }
}

export const binom = {
    // This implementation comes from:
    // http://home.ubalt.edu/ntsbarsh/business-stat/otherapplets/ConfIntPro.htm
    test: function (m, n, level) {
        if (typeof(level)=='undefined') {
            level = 0.95;
        }

        var P1 = Math.round((m / n) * 100000) / 100000;
        var P = m / n;
        var nu1 = 2 * (n - m + 1);
        var nu2 = 2 * m;
        var AL = (1 - level) / 2;
        var F1 = this.AFishF(AL, nu1, nu2);
        var NUM1 = m;
        var DEN1 = (n - m + 1) * F1 + m;
        var LL1 = (NUM1 / DEN1);
        var nup1 = nu2 + 2;
        var nup2 = nu1 - 2;
        var F2 = this.AFishF(AL, nup1, nup2)
        var NUM2 = (m + 1) * F2;
        var DEN2 = (n - m) + (m + 1) * F2;
        var UL1 = (NUM2 / DEN2);

        var lower = Math.round(100000 * LL1) / 100000;
        var upper = Math.round(100000 * UL1) / 100000;

        return {
            estimate: P1,
            lower: lower,
            upper: upper
        };

    },

    FishF: function (f, n1, n2) {
        var Pi = Math.PI;
        var PiD2 = Pi / 2;
        var PiD4 = Pi / 4;
        var Pi2 = 2 * Pi;
        var e = 2.718281828459045235;
        var e10 = 1.105170918075647625;
        var Deg = 180 / Pi;
        var x = n2 / (n1 * f + n2);
        if ((n1 % 2) == 0) { 
            return this.StatCom(1 - x, n2, n1 + n2 - 4, n2 - 2) * Math.pow(x, n2 / 2);
        }
        if ((n2 % 2) == 0) { 
            return 1 - this.StatCom(x, n1, n1 + n2 - 4, n1 - 2) * Math.pow(1 - x, n1 / 2);
        }
        var th = Math.atan(Math.sqrt(n1 * f / n2)); 
        var a = th / PiD2; var sth = Math.sin(th);
        var cth = Math.cos(th);
        if (n2 > 1) { 
            a = a + sth * cth * this.StatCom(cth * cth, 2, n2 - 3, -1) / PiD2;
        }
        if (n1 == 1) { 
            return 1 - a;
        }
        var c = 4 * this.StatCom(sth * sth, n2 + 1, n1 + n2 - 4, n2 - 2) * sth * Math.pow(cth, n2) / Pi;
        if (n2 == 1) { 
            return 1 - a + c / 2; 
        }
        var k = 2; 
        while (k <= (n2 - 1) / 2) {
            c = c * k / (k - .5); 
            k = k + 1;
        }

        return 1 - a + c;
    },

    StatCom: function (q, i, j, b) {
        var zz = 1; 
        var z = 1; 
        var k = i; 
        while (k <= j) { 
            zz = zz * q * k / (k - b); 
            z = z + zz; 
            k = k + 2; 
        }
        return z;
    },

    AFishF: function (p, n1, n2) {
        var v = 0.5; 
        var dv = 0.5; 
        var f = 0;

        while (dv > 1e-10) { 
            f = 1 / v - 1; 
            dv = dv / 2; 

            if (this.FishF(f, n1, n2) > p) { 
                v = v - dv;
            } else { 
                v = v + dv;
            } 
        }
        return f;
    }
};
