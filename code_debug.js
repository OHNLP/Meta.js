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
        {treat1: 'A', treat2: 'B'}
    ]);

    console.log(ret);

    // test 2
    ret = metajs.calc_n_graphs([
        {treat1: 'A', treat2: 'B'},
        {treat1: 'A', treat2: 'C'},
        {treat1: 'A', treat2: 'D'},
        {treat1: 'A', treat2: 'E'},
    ]);

    console.log(ret);

    // test 3
    ret = metajs.calc_n_graphs([
        {treat1: 'A', treat2: 'B'},
        {treat1: 'K', treat2: 'C'},
        {treat1: 'A', treat2: 'D'},
        {treat1: 'A', treat2: 'E'},
        {treat1: 'F', treat2: 'G'},
    ]);

    // test 5
    ret = metajs.calc_n_graphs([
        {treat1: 'A', treat2: 'B'},
        {treat1: 'K', treat2: 'C'},
        {treat1: 'H', treat2: 'D'},
        {treat1: 'I', treat2: 'E'},
        {treat1: 'F', treat2: 'G'},
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

function debug_netmeta() {
    var rs = [
        {study: 'TRA', sm: 0.62, lower: 0.46, upper: 0.82, treat1: 'CaboNivo', treat2: 'Suni', year: 2020},
        {study: 'TRB', sm: 0.54, lower: 0.46, upper: 0.63, treat1: 'NivoIpi', treat2: 'Suni', year: 2021},
        {study: 'TRC', sm: 0.68, lower: 0.58, upper: 0.81, treat1: 'AteBev', treat2: 'Suni', year: 2022},
        {study: 'TRD', sm: 0.80, lower: 0.62, upper: 1.03, treat1: 'Pazo', treat2: 'Suni', year: 2022},
        {study: 'TRE', sm: 0.74, lower: 0.45, upper: 1.20, treat1: 'AteBev', treat2: 'Suni', year: 2022},
        {study: 'TRF', sm: 1.12, lower: 0.91, upper: 1.38, treat1: 'PemAxi', treat2: 'Suni', year: 2022},
    ];

    var ret = metajs.netmeta(rs, {});

    console.log(ret);
};

debug_netmeta();
