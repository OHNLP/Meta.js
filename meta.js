/**
 * Meta-analysis in JavaScript
 * 
 * Please make sure the num.js and math.js are loaded
 */
var metajs = {
    backtransf: function(x, sm, value, n) {
        return x;
    },

    asin2p: function(x, n, value) {
        if (typeof(n) == 'undefined') {
            n = null;
        }
        if (typeof(value) == 'undefined') {
            value = 'mean';
        }

        var minimum = Math.asin(0);
        var maximum = Math.asin(1);

        if (n != null) {
            minimum = 0.5 * (Math.asin(Math.sqrt(0 / (n + 1))) + Math.asin(Math.sqrt((0 + 1) / (n + 1))))
            maximum = 0.5 * (Math.asin(Math.sqrt(n / (n + 1))) + Math.asin(Math.sqrt((n + 1) / (n + 1))))
        }
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

        if (!params.hasOwnProperty('input_format')) {
            params['input_format'] = 'PRIM_CAT_RAW';
        }

        if (!params.hasOwnProperty('sm')) {
            params['sm'] = 'PLOGIT';
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
            n: []
        };
        // a mapping from ds index to rs index
        var d2r = {};
        for (let i = 0; i < rs.length; i++) {
            const r = rs[i];
            // check the r
            if (r[0] > r[1]) {
                // it's not possible that event > n
                continue;

            } else if (r[0] == 0) {
                // zero event??? increase both
                ds.e.push( incr );
                ds.n.push( r[1] );

            } else {
                // for most case
                ds.e.push( r[0] );
                ds.n.push( r[1] );
            }
            // add this mapping
            d2r[ds.e.length - 1] = i;
        }
        // double check the length of records
        if (ds.e.length == 0) {
            // what???
        } else if (ds.e.length == 1) {
            // what??? only one study?
        } else {
            // ok, more than 1
        }

        // convert the e and n to numjs format
        ds.e = nj.array(ds.e);
        ds.n = nj.array(ds.n);
        // one for calcualtion
        ONE = nj.ones(ds.e.size);

        ///////////////////////////////////////////////////
        // (6) Subset analysis
        ///////////////////////////////////////////////////


        ///////////////////////////////////////////////////
        // (7) Calculate results for each study
        ///////////////////////////////////////////////////
        var TE = null;
        var seTE = null;

        if (params.sm == 'PLOGIT') {
            TE = nj.log(
                ds.e.divide(
                    ds.n.subtract(ds.e)
                )
            );
            seTE = (ONE.divide(ds.e).add(
                ONE.divide(ds.n.subtract(ds.e))
            )).pow(0.5);

        } else if (params.sm == 'PFT') {
            TE = math.dotMultiply(
                0.5,
                math.add(
                    math.asin(math.sqrt(math.dotDivide(e, math.add(n, 1)))),
                    math.asin(math.sqrt(math.dotDivide(math.add(e, 1), math.add(n, 1))))
                )
            );
            seTE = math.sqrt(
                math.dotDivide(
                    1,
                    math.add(
                        2,
                        math.dotMultiply(
                            4,
                            n
                        )
                    )
                )
            );
        }

        var SM = [];
        var SM_lower = [];
        var SM_upper = [];

        // var SM_RS = ds.e.tolist().map((e,i)=>binom.test(e, ds.n.get(i)));
        for (let i = 0; i < ds.e.size; i++) {
            var _e = ds.e.get(i);
            var _n = ds.n.get(i);
            var r = binom.test(_e, _n);

            SM.push(r.estimate);
            SM_lower.push(r.lower);
            SM_upper.push(r.upper);
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
        var w_fixed = ONE.divide(seTE.pow(2));
        var wp_fixed = w_fixed.divide(w_fixed.sum());

        var TE_fixed = TE.multiply(wp_fixed).sum();
        var SM_fixed = this.expit(TE_fixed);
        var seTE_fixed = math.sqrt( 1 / w_fixed.sum());

        var SM_fixed_lower = this.expit(TE_fixed - 1.96 * seTE_fixed);
        var SM_fixed_upper = this.expit(TE_fixed + 1.96 * seTE_fixed);

        fixed = {
            TE: TE_fixed,
            seTE: seTE_fixed,
            w: w_fixed.tolist(),
            wp: wp_fixed.tolist(),

            SM: SM_fixed,
            SM_lower: SM_fixed_lower,
            SM_upper: SM_fixed_upper
        }
        ///////////////////////////////
        // Heteroginity
        ///////////////////////////////
        var heterogeneity = this.heterogeneity_by_DL(TE, seTE);

        ///////////////////////////////
        // Random effect model
        ///////////////////////////////

        // TODO

        ///////////////////////////////////////////////////
        // (9) Finalize return object
        ///////////////////////////////////////////////////
        var ret = {
            ds: {
                e: ds.e.tolist(),
                n: ds.n.tolist()
            },
            heterogeneity: heterogeneity,

            TE: TE.tolist(),
            seTE: seTE.tolist(),
            SM: SM,
            SM_lower: SM_lower,
            SM_upper: SM_upper,

            // the MA result
            fixed: fixed,
            random: random,

            // the settings
            params: params
        }

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
            params['method'] = 'Inverse';
        }

        if (!params.hasOwnProperty('method_ci')) {
            params['method_ci'] = 'CP';
        }

        if (!params.hasOwnProperty('method_tau')) {
            params['method_tau'] = 'DL';
        }

        ///////////////////////////////////////////////////
        // (2) Read data
        ///////////////////////////////////////////////////


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
            Nc: []
        };
        // a mapping from ds index to rs index
        var d2r = {};
        for (let i = 0; i < rs.length; i++) {
            const r = rs[i];
            // check the r
            if (r[0] > r[1]) {
                // it's not possible that event > n
                continue;

            } else if (r[0] == 0 || r[2] == 0) {
                // zero event??? increase both
                ds.Et.push( r[0] + incr );
                ds.Nt.push( r[1] + incr );
                ds.Ec.push( r[2] + incr );
                ds.Nc.push( r[3] + incr );

            } else {
                // for most case
                ds.Et.push( r[0] );
                ds.Nt.push( r[1] );
                ds.Ec.push( r[2] );
                ds.Nc.push( r[3] );
            }
            // add this mapping
            d2r[ds.Et.length - 1] = i;
        }
        // double check the length of records
        if (ds.Et.length == 0) {
            // what???
        } else if (ds.Et.length == 1) {
            // what??? only one study?
        } else {
            // ok, more than 1
        }

        // convert the e and n to numjs format
        ds.Et = nj.array(ds.Et);
        ds.Nt = nj.array(ds.Nt);
        ds.Ec = nj.array(ds.Ec);
        ds.Nc = nj.array(ds.Nc);

        // for the calculation
        ds.n11 = ds.Et;
        ds.n21 = ds.Ec;
        ds.n1_ = ds.Nt;
        ds.n2_ = ds.Nc;
        ds.n__ = ds.Nt.add(ds.Nc);
        ds.n12 = ds.Nt.subtract(ds.Et);
        ds.n22 = ds.Nc.subtract(ds.Ec);
        ds.n_1 = ds.n11.add(ds.n21);
        ds.n_2 = ds.n12.add(ds.n22);

        // one for calcualtion
        ONE = nj.ones(ds.Et.size);

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
                TE = nj.log(
                    ds.n11.multiply(ds.n22).divide(
                        ds.n12.multiply(ds.n21)
                    )
                );
                seTE = (ONE.divide(ds.n11)
                    .add(ONE.divide(ds.n12))
                    .add(ONE.divide(ds.n21))
                    .add(ONE.divide(ds.n22))
                ).pow(0.5);

                // backtransf
                SM = nj.exp(TE);
                SM_lower = nj.exp(TE.subtract(seTE.multiply(1.96)));
                SM_upper = nj.exp(TE.add(seTE.multiply(1.96)));
            }

        } else if (params.sm == 'RR') {
            TE = nj.log(
                (ds.n11.divide(ds.n1_)).divide(
                    ds.n21.divide(ds.n2_)
                )
            );
            // Hartung & Knapp (2001), Stat Med, equation (18)
            seTE = (
                ONE.divide(ds.n11)
                .subtract(ONE.divide(ds.n1_))
                .add(ONE.divide(ds.n21))
                .subtract(ONE.divide(ds.n2_))
            ).pow(0.5);

            // backtransf 
            SM = nj.exp(TE);
            SM_lower = nj.exp(TE.subtract(seTE.multiply(1.96)));
            SM_upper = nj.exp(TE.add(seTE.multiply(1.96)));
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
        if (['MH', 'Inverse'].includes(params.method)) {
            if (params.sm == 'OR') {
                var A = ds.n11.multiply(ds.n22).divide(ds.n__);
                var B = ds.n11.add(ds.n22).divide(ds.n__);
                var C = ds.n12.multiply(ds.n21).divide(ds.n__);
                var D = ds.n12.add(ds.n21).divide(ds.n__);

                var w_fixed = C;
                var wp_fixed = w_fixed.divide(w_fixed.sum());

                var TE_fixed = math.log(A.sum() / C.sum());
                var seTE_fixed = math.sqrt(
                    (1 / (2 * A.sum()**2)) * (
                        A.multiply(B).sum() 
                        +
                        math.exp(TE_fixed) * (
                            B.multiply(C).sum() +
                            A.multiply(D).sum()
                        ) 
                        +
                        math.exp(TE_fixed)**2 * C.multiply(D).sum()
                    )
                );

                var SM_fixed = math.exp(TE_fixed);
                var SM_fixed_lower = math.exp(TE_fixed - 1.96 * seTE_fixed);
                var SM_fixed_upper = math.exp(TE_fixed + 1.96 * seTE_fixed);

                fixed = {
                    TE: TE_fixed,
                    seTE: seTE_fixed,
                    w: w_fixed.tolist(),
                    wp: wp_fixed.tolist(),

                    SM: SM_fixed,
                    SM_lower: SM_fixed_lower,
                    SM_upper: SM_fixed_upper
                }

            } else if (params.sm == 'RR') {
                var D = (ds.n1_.multiply(ds.n2_).multiply(ds.n_1)).subtract(
                    ds.n11.multiply(ds.n21).multiply(ds.n__)
                ).divide(ds.n__.pow(2));
                var R = ds.n11.multiply(ds.n2_).divide(ds.n__);
                var S = ds.n21.multiply(ds.n1_).divide(ds.n__);

                var w_fixed = S;
                var wp_fixed = w_fixed.divide(w_fixed.sum());

                var TE_fixed = math.log(R.sum() / S.sum());
                var seTE_fixed = math.sqrt(D.sum() / (R.sum() * S.sum()));

                var SM_fixed = math.exp(TE_fixed);
                var SM_fixed_lower = math.exp(TE_fixed - 1.96 * seTE_fixed);
                var SM_fixed_upper = math.exp(TE_fixed + 1.96 * seTE_fixed);

                fixed = {
                    TE: TE_fixed,
                    seTE: seTE_fixed,
                    w: w_fixed.tolist(),
                    wp: wp_fixed.tolist(),

                    SM: SM_fixed,
                    SM_lower: SM_fixed_lower,
                    SM_upper: SM_fixed_upper
                }
            }
        }
        
        ///////////////////////////////
        // Heteroginity
        ///////////////////////////////
        var heterogeneity = this.heterogeneity_by_DL(TE, seTE);

        ///////////////////////////////
        // Random effect model
        ///////////////////////////////

        // TODO

        ///////////////////////////////////////////////////
        // (9) Finalize return object
        ///////////////////////////////////////////////////
        var ret = {
            ds: {
                Et: ds.Et.tolist(),
                Nt: ds.Nt.tolist(),
                Ec: ds.Ec.tolist(),
                Nc: ds.Nc.tolist()
            },
            heterogeneity: heterogeneity,

            // each study
            TE: TE.tolist(),
            seTE: seTE.tolist(),
            SM: SM.tolist(),
            SM_lower: SM_lower.tolist(),
            SM_upper: SM_upper.tolist(),

            // the MA result
            fixed: fixed,
            random: random,

            // the settings
            params: params
        }

        return ret;
    },


    expit: function(v) {
        if (typeof(v) == 'number') {
            return this.expit_value(v);

        } else if (v.hasOwnProperty('_size')) {
            // which means it is a mathjs matrix
            return this.expit_mathjs(v);

        } else if (v.hasOwnProperty('size')) {
            // which means it is a nj obj
            return this.expit_nj(v);

        } else if (v.hasOwnProperty('length')) {
            // which ma
            return this.expit_mathjs(v);

        } else {
            return null;
        }
    },

    expit_value: function(v) {
        return 0.5 + 0.5 * math.tanh(v*0.5);
    },

    expit_mathjs: function(v) {
        return math.add(
            math.multiply(
                0.5,
                math.tanh(
                    math.multiply(
                        v,
                        0.5
                    )
                )
            ),
            0.5
        );
    },

    expit_nj: function(V) {
        return nj.tanh(V.multiply(0.5)).add(1).multiply(0.5);
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
            TE: Math.exp(tres.TE),
            lower: Math.max(Math.exp(tres.lower), 1),
            upper: Math.max(Math.exp(tres.upper), 1)
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

        var statistic = math.divide(TE, seTE);

        var pval = 2 * pnorm.compute(math.abs(statistic), 0, 1, false);

        return {
            TE: TE,
            seTE: seTE,
            pval: pval,
            level: level,
            lower: lower,
            upper: upper,
            statistic: statistic
        };
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
    heterogeneity_by_DL: function(TE, seTE) {
        var k = TE.size;
        var p = 1;
        var ONE = nj.ones(TE.size);

        var vi = seTE.pow(2);
        var wi = ONE.divide(vi);
        var W = nj.array(
            math.diag(
                wi.tolist()
            )
        ).tolist();
        var X = math.ones(TE.size, 1);
        var Y = math.transpose([TE.tolist()]);
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

        return {
            I2: I2.TE,
            tau2: tau2,
            se_tau2: se_tau2,
            pval_Q: pval_Q
        };
    }

}

var binom = {
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
            c = c * k / (k - .5); k = k + 1 
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


var pchisq = {
    // https://www.math.ucla.edu/~tom/distributions/chisq.html
    log_gamma: function(Z) {
        with (Math) {
            var S=1+76.18009173/Z-86.50532033/(Z+1)+24.01409822/(Z+2)-1.231739516/(Z+3)+.00120858003/(Z+4)-.00000536382/(Z+5);
            var LG= (Z-.5)*log(Z+4.5)-(Z+4.5)+log(S*2.50662827465);
        }
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
        while (Math.abs((A1-AOLD)/A1)>.00001) {
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
        var Prob=Math.exp(A*Math.log(X)-X- this.log_gamma(A))*A1;
    
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
        G = G*Math.exp(A*Math.log(X)-X- this.log_gamma(A));

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
        
        p = this.gamma_cdf(q/2, df/2);
        p = Math.round(p*100000) / 100000;

        if (lower_tail) {

        } else {
            p = 1 - p;
            p = Math.round(p*100000) / 100000;
        }
        
        return p;
    }
}


var pnorm = {
    normal_cdf: function(X){
        var T=1/(1+.2316419*Math.abs(X));
        var D=.3989423*Math.exp(-X*X/2);
        var p=D*T*(.3193815+T*(-.3565638+T*(1.781478+T*(-1.821256+T*1.330274))));
        if (X>0) {
            p=1-p
        }
        return p
    },

    compute: function(q, mean, sd, lower_tail) {
        if (typeof(mean) == 'undefined') {
            mean = 0;
        }
        if (typeof(sd) == 'undefined') {
            sd = 1;
        }
        if (typeof(lower_tail) == 'undefined') {
            lower_tail = true;
        }
        if (sd<0) {
            return NaN;
        }
        if (sd==0) {
            if (q < mean){
                p = 0
            } else {
                p = 1
            }

        } else {
            p = this.normal_cdf((q-mean)/sd);
            p = Math.round(100000*p)/100000;
        }

        if (lower_tail) {

        } else {
            p = 1 - p;
            p = Math.round(100000*p)/100000;
        }
        
        return p;
    }
}
