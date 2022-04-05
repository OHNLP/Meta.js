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

debug_metabin();
// debug_metaprop();
// debug_csv_read();