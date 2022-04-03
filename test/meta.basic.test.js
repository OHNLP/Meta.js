'use strict';

import { metajs } from '../src/meta.js';
import assert from 'assert'

describe('testing metajs base functions', () => {
    it('backtransf 1 OR null null should be 1', () => {

        assert.equal(
            metajs.backtransf(1, 'OR', null, null), 
            1
        );
    });

    // related to ci95
    it('ci95(0.5, 0.2) should be ', () => {
        assert.deepEqual(
            metajs.ci95(0.5, 0.2), 
            {
                "TE": 0.5,
                "seTE": 0.2,
                "pval": 0.01242,
                "level": 0.95,
                "lower": 0.108,
                "upper": 0.892,
                "statistic": 2.5
            } 
        );
    });

    // just an empty test
    it('at least one pass', () => {
        assert.equal(
            1,
            1
        );
    });
});


describe('testing metajs.metabin functions', () => {
    it('simple records, fixed OR by MH should be 1.697488', () => {
        var d = {
            rs: [
                [12,393,2, 396, 'S1', 'T','C'],
                [24,230,24,281, 'S2', 'T','C'],
            ],
            vs: {
                fixed: {
                    SM: 1.697488,
                    SM_lower: 0.999475,
                    SM_upper: 2.882978
                }
            }
        };

        var vals = metajs.metabin(d.rs);
        assert.equal(
            metajs.tfxd6(vals.fixed.SM),
            d.vs.fixed.SM
        );
        assert.equal(
            metajs.tfxd6(vals.fixed.SM_lower),
            d.vs.fixed.SM_lower
        );
        assert.equal(
            metajs.tfxd6(vals.fixed.SM_upper),
            d.vs.fixed.SM_upper
        );
    });

    it('simple with 0 event, fixed OR by MH should be 1.13', () => {
        var d = {
            rs: [
                [0 ,393,2, 396, 'S1', 'T','C'],
                [24,230,24,281, 'S2', 'T','C'],
            ],
            vs: {
                fixed: ['1.13', '0.64', '2.00']
            }
        };

        var vals = metajs.metabin(d.rs);
        assert.deepEqual(
            [
                vals.fixed.SM.toFixed(2), 
                vals.fixed.SM_lower.toFixed(2), 
                vals.fixed.SM_upper.toFixed(2)
            ],
            d.vs.fixed
        );
    });
});