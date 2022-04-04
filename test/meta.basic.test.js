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


describe('testing metajs.metaprop functions', () => {
    it('simple case, fixed incd by Inverse should be 0.15', () => {
        var d = {
            rs: [
                [2, 20,  'S1'],
                [5, 90,  'S2'],
                [20,100, 'S3'],
            ],
            vs: {
                fixed: ['0.15', '0.10', '0.21'],
                random: ['0.11', '0.04', '0.26']
            }
        }

        var vals = metajs.metaprop(d.rs, { sm: 'PLOGIT' });
        assert.deepEqual(
            [
                vals.fixed.SM.toFixed(2), 
                vals.fixed.SM_lower.toFixed(2), 
                vals.fixed.SM_upper.toFixed(2)
            ],
            d.vs.fixed
        );
    });

    it('simple case, random incd by Inverse should be 0.11', () => {
        var d = {
            rs: [
                [2, 20,  'S1'],
                [5, 90,  'S2'],
                [20,100, 'S3'],
            ],
            vs: {
                fixed: ['0.15', '0.10', '0.21'],
                random: ['0.11', '0.04', '0.26']
            }
        }

        var vals = metajs.metaprop(d.rs, { sm: 'PLOGIT' });
        assert.deepEqual(
            [
                vals.random.SM.toFixed(2), 
                vals.random.SM_lower.toFixed(2), 
                vals.random.SM_upper.toFixed(2)
            ],
            d.vs.random
        );
    })
});


describe('testing metajs.metabin functions', () => {
    it('simple case, fixed OR by MH should be 1.697488', () => {
        var d = {
            rs: [
                [12,393,2, 396, 'S1', 'T','C'],
                [24,230,24,281, 'S2', 'T','C'],
            ],
            vs: {
                fixed: ['1.697', '0.999', '2.883']
            }
        };

        var vals = metajs.metabin(d.rs);
        assert.deepEqual(
            [
                vals.fixed.SM.toFixed(3), 
                vals.fixed.SM_lower.toFixed(3), 
                vals.fixed.SM_upper.toFixed(3)
            ],
            d.vs.fixed
        );
    });

    it('simple case, fixed RR by MH should be 1.63', () => {
        var d = {
            rs: [
                [12,393,2, 396, 'S1', 'T','C'],
                [24,230,24,281, 'S2', 'T','C'],
            ],
            vs: {
                fixed: ['1.63', '1.00', '2.66']
            }
        };

        var vals = metajs.metabin(d.rs, { sm: 'RR' });
        assert.deepEqual(
            [
                vals.fixed.SM.toFixed(2), 
                vals.fixed.SM_lower.toFixed(2), 
                vals.fixed.SM_upper.toFixed(2)
            ],
            d.vs.fixed
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

    it('simple with 0 event, fixed RR by MH should be 1.12', () => {
        var d = {
            rs: [
                [0 ,393,2, 396, 'S1', 'T','C'],
                [24,230,24,281, 'S2', 'T','C'],
            ],
            vs: {
                fixed: ['1.12', '0.66', '1.88']
            }
        };

        var vals = metajs.metabin(d.rs, { sm: 'RR' });
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