'use strict';

import {create, all, random} from 'mathjs'
const config = { }
const math = create(all, config);

import { metajs } from './src/meta.js';
import Papa from 'papaparse';
import fs from 'fs';
import * as dfd from "danfojs-node";

function debug_csv_read() {

    // locate the input stream
    let f = fs.createReadStream('./testsets/test_input.csv');
    let fo = fs.createReadStream('./testsets/test_result_OR.csv');

    var rs1 = null;
    var rs2 = null;

    Papa.parse(f, {
        header: true,
        dynamicTyping: true,
        comments: "#",
        complete: function (rs) {
            console.log('* loaded ' + rs.data.length + ' records in input testset');
            rs1 = rs;
            Papa.parse(fo, {
                header: true,
                dynamicTyping: true,
                comments: "#",
                complete: function (rs) {
                    console.log('* loaded ' + rs.data.length + ' records in result testset');
                    rs2 = rs;

                    // ok, now let's test
                    do_something()
                }
            });
        }
    });

    function do_something() {
        // create a dataframe
        var df = new dfd.DataFrame(rs1.data);
        var oc_names = df['outcome'].unique().values;

        // get the dataframe for results
        var dfr = new dfd.DataFrame(rs2.data);

        function _tfx(v) {
            if (isNaN(v)) {
                return 'NA';
            }
            if (v == null) {
                return 'NA'
            }
            if (v == 'NA') {
                return 'NA';
            }
            return v.toFixed(2);
        }

        function pad(pad, str, padLeft) {
            if (typeof str === 'undefined')
                return pad;
            if (padLeft) {
                return (pad + str).slice(-pad.length);
            } else {
                return (str + pad).substring(0, pad.length);
            }
        }

        console.log('*--------------------------------');
        console.log('* fixed results:');
        console.log('*           outcome           \tmeta.js\t\t\t| R.meta\t\t|')
        console.log('*--------------------------------');
        for (let i = 0; i < oc_names.length; i++) {
            const ocn = oc_names[i];
            
            var vals = df.loc({
                rows: df['outcome'].eq(ocn), 
                columns: ['Et', 'Nt', 'Ec', 'Nc', 'study']
            }).values;
    
            var gtrs = dfr.loc({
                rows: dfr['outcome'].eq(ocn), 
                columns: ['TE.fixed', 'lower.fixed', 'upper.fixed']
            }).values[0];
    
            var rst = metajs.metabin(vals, {
                sm: 'OR'
            });

            console.log(
                "* " + pad('                      ', ocn, true) + '\t' + 
                _tfx(rst.fixed.TE) + '\t' + 
                _tfx(rst.fixed.TE_lower) + '\t' + 
                _tfx(rst.fixed.TE_upper) + '\t| ' +

                _tfx(gtrs[0]) + '\t' + 
                _tfx(gtrs[1]) + '\t' + 
                _tfx(gtrs[2]) + '\t|' 

            );
        }


        console.log('*--------------------------------');
        console.log('* random results:');
        console.log('*           outcome           \tmeta.js\t\t\t| R.meta\t\t|')
        console.log('*--------------------------------');
        for (let i = 0; i < oc_names.length; i++) {
            const ocn = oc_names[i];
            
            var vals = df.loc({
                rows: df['outcome'].eq(ocn), 
                columns: ['Et', 'Nt', 'Ec', 'Nc', 'study']
            }).values;
    
            var gtrs = dfr.loc({
                rows: dfr['outcome'].eq(ocn), 
                columns: ['TE.random', 'lower.random', 'upper.random']
            }).values[0];
    
            var rst = metajs.metabin(vals, {
                sm: 'OR'
            });

            console.log(
                "* " + pad('                      ', ocn, true) + '\t' + 
                _tfx(rst.random.TE) + '\t' + 
                _tfx(rst.random.TE_lower) + '\t' + 
                _tfx(rst.random.TE_upper) + '\t| ' +

                _tfx(gtrs[0]) + '\t' + 
                _tfx(gtrs[1]) + '\t' + 
                _tfx(gtrs[2]) + '\t|' 

            );
        }
    }
}

function debug_metabin() {
    var rs = [
        [12,393,2,396, 'S1', 'T','C'],
        [24,230,24,281, 'S2', 'T','C'],
    ]
    
    var rst = metajs.metabin(
        rs,
        {
            'sm': 'OR'
        }
    );
    
    console.log(rst);
}

function debug_metaprop() {
    var rs = [
        [2, 20,  'S1'],
        [5, 90,  'S2'],
        [20,100, 'S3'],
    ];

    var rst = metajs.metaprop(
        rs, 
        {
            'sm': 'PFT'
        }
    );

    console.log(rst);
}

function debug_calc_n_graphs() {
    var ret = metajs.calc_n_graphs([
        {t1: 'A', t2: 'B'}
    ]);

    console.log(ret);

    // test 2
    ret = metajs.calc_n_graphs([
        {t1: 'A', t2: 'B'},
        {t1: 'A', t2: 'C'},
        {t1: 'A', t2: 'D'},
        {t1: 'A', t2: 'E'},
    ]);

    console.log(ret);

    // test 3
    ret = metajs.calc_n_graphs([
        {t1: 'A', t2: 'B'},
        {t1: 'K', t2: 'C'},
        {t1: 'A', t2: 'D'},
        {t1: 'A', t2: 'E'},
        {t1: 'F', t2: 'G'},
    ]);

    // test 5
    ret = metajs.calc_n_graphs([
        {t1: 'A', t2: 'B'},
        {t1: 'K', t2: 'C'},
        {t1: 'H', t2: 'D'},
        {t1: 'I', t2: 'E'},
        {t1: 'F', t2: 'G'},
    ]);

    console.log(ret);
}

// sample data
// CaboNivo	Suni	0.62	0.46	0.82
// NivoIpi	Suni	0.54	0.46	0.63
// AteBev	Suni	0.68	0.58	0.81
// Pazo	    Suni	0.80	0.62	1.03
// AteBev	Suni	0.74	0.45	1.2
// PemAxi	Suni	1.12	0.91	1.38

function _debug_netmeta_test1() {
    
    // test 1
    var rs = [
        {study: 'SA', sm: 0.62, lower: 0.46, upper: 0.82, t1: 'CaboNivo', t2: 'Suni', year: 2020},
        {study: 'SB', sm: 0.54, lower: 0.46, upper: 0.63, t1: 'NivoIpi', t2: 'Suni', year: 2021},
        {study: 'SC', sm: 0.68, lower: 0.58, upper: 0.81, t1: 'AteBev', t2: 'Suni', year: 2022},
        {study: 'SD', sm: 0.80, lower: 0.62, upper: 1.03, t1: 'Pazo', t2: 'Suni', year: 2022},
        {study: 'SE', sm: 0.74, lower: 0.45, upper: 1.20, t1: 'AteBev', t2: 'Suni', year: 2022},
        {study: 'SF', sm: 1.12, lower: 0.91, upper: 1.38, t1: 'PemAxi', t2: 'Suni', year: 2022},
    ];

    var nma = metajs.netmeta(rs, {});

    // the league table
    metajs.print_league_table(nma);

    // the rank list
    var rank = metajs.netrank(nma);
    metajs.print_network_rank(nma, rank);

};

function _debug_netmeta_test2() {
    // test 2
    var rs = [
        {t1: 'E_ADT', t2: 'ADT', sm: 0.72, lower: 0.47, upper: 1.09, study: 'ENZAMET' },
        {t1: 'APA_ADT', t2: 'ADT', sm: 0.4, lower: 0.15, upper: 1.03, study: 'TITAN' },
        {t1: 'D_ADT', t2: 'ADT', sm: 0.83, lower: 0.47, upper: 1.47, study: 'GETUG_AFU15' },
        {t1: 'DARO_D_ADT', t2: 'D_ADT', sm: 0.605, lower: 0.348, upper: 1.052, study: 'ARASENS' },
    ];

    // the league table
    var nma = metajs.netmeta(rs, {});
    metajs.print_league_table(nma);

    // 
    var rank = metajs.netrank(nma);
    metajs.print_network_rank(nma, rank);
}

function debug_netmeta() {
    _debug_netmeta_test1();
}

debug_netmeta();
