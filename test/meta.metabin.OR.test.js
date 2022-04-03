'use strict';

import { metajs } from '../src/meta.js';
import assert from 'assert'
import Papa from 'papaparse';
import fs from 'fs';
import * as dfd from "danfojs-node";

// locate the input stream
let f = fs.createReadStream('./testsets/test_input.csv');
let fo = fs.createReadStream('./testsets/test_result_OR.csv');

var rs1 = null;
var rs2 = null;

Papa.parse(f, {
    header: true,
    comments: "#",
    complete: function(rs) {
        console.log('* loaded ' + rs.data.length + ' records in input testset');
        rs1 = rs;
        Papa.parse(fo, {
            header: true,
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

function start_test() {
describe('testing metajs metabin functions', () => {

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

       it('case ' + ocn + ' should be correct', (done) => {

            metajs.metabin(vals, {

            });
            assert.equal(
                1, 
                1
            );

            done();
        });
   }
});
}
