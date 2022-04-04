'use strict';

import { metajs } from '../src/meta.js';
import assert from 'assert'
import Papa from 'papaparse';
import fs from 'fs';
import * as dfd from "danfojs-node";

// locate the input stream
let f = fs.createReadStream('./testsets/test_input.csv');
let fo = fs.createReadStream('./testsets/test_result_PFT.csv');

var rs1 = null;
var rs2 = null;

Papa.parse(f, {
    header: true,
    dynamicTyping: true,
    comments: "#",
    complete: function(rs) {
        console.log('* loaded ' + rs.data.length + ' records in input testset');
        rs1 = rs;
        Papa.parse(fo, {
            header: true,
            dynamicTyping: true,
            comments: "#",
            complete: function(rs) {
                console.log('* loaded ' + rs.data.length + ' records in result testset');
                rs2 = rs;

                // ok, now let's test
                start_test();
            }
        });
    }
});

function _tf2(v) {
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

function start_test() {
describe('testing metajs metaprop for PFT functions', () => {

   // create a dataframe
   var df = new dfd.DataFrame(rs1.data);
   var oc_names = df['outcome'].unique().values;

   // get the dataframe for results
   var dfr = new dfd.DataFrame(rs2.data);

   for (let i = 0; i < oc_names.length; i++) {
       const ocn = oc_names[i];
       
       var vals = df.loc({
           rows: df['outcome'].eq(ocn), 
           columns: ['Et', 'Nt', 'Ec', 'Nc', 'study']
       }).values;

       var gtrs = dfr.loc({
           rows: dfr['outcome'].eq(ocn), 
           columns: ['TE.fixed', 'lower.fixed', 'upper.fixed',
                     'TE.random', 'lower.random', 'upper.random']
        }).values[0];

       it('case ' + ocn + ', fixed PFT should be ' + _tf2(gtrs[0]) + '(' +_tf2(gtrs[1]) + ',' + _tf2(gtrs[2]) + ')', (done) => {
            var rst = metajs.metaprop(vals, {
                sm: 'PFT'
            });

            assert.deepEqual(
                [_tf2(rst.fixed.TE), _tf2(rst.fixed.TE_lower), _tf2(rst.fixed.TE_upper)], 
                [_tf2(gtrs[0]), _tf2(gtrs[1]), _tf2(gtrs[2])]
            );

            done();
        });

        it('case ' + ocn + ', random PFT should be ' + _tf2(gtrs[3]) + '(' +_tf2(gtrs[4]) + ',' + _tf2(gtrs[5]) + ')', (done) => {
            var rst = metajs.metaprop(vals, {
                sm: 'PFT'
            });

            assert.deepEqual(
                [_tf2(rst.random.TE), _tf2(rst.random.TE_lower), _tf2(rst.random.TE_upper)], 
                [_tf2(gtrs[3]), _tf2(gtrs[4]), _tf2(gtrs[5])]
            );

            done();
        });
   }
});
}
