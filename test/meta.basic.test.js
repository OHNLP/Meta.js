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


describe('testing metajs network functions', () => {
    it('calc_n_graphs one graph should be 1', () => {
        assert.equal(
            metajs.calc_n_graphs([
                {treat1: 'A', treat2: 'B'}
            ]).length, 
            1
        );
    });

    it('calc_n_graphs two graphs should be 2', () => {
        assert.equal(
            metajs.calc_n_graphs([
                {treat1: 'A', treat2: 'B'},
                {treat1: 'C', treat2: 'D'},
            ]).length, 
            2
        );
    });

    it('calc_n_graphs a lot of edges should be 1', () => {
        assert.equal(
            metajs.calc_n_graphs([
                {treat1: 'Nivo', treat2: 'Suni'},
                {treat1: 'Cabo', treat2: 'Suni'},
                {treat1: 'CaboNivo', treat2: 'Suni'},
                {treat1: 'Treat', treat2: 'Suni'},
            ]).length, 
            1
        );
    });
});


describe('testing metajs.metaprop functions', () => {
    it('simple case, fixed incd/PLOGIT by Inverse should be 0.15 (0.10, 0.21)', () => {
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

    it('simple case, random incd/PLOGIT by Inverse should be 0.11 (0.04, 0.26)', () => {
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
    });

    it('simple case, fixed incd/PFT by Inverse should be 0.12 (0.07, 0.17)', () => {
        var d = {
            rs: [
                [2, 20,  'S1'],
                [5, 90,  'S2'],
                [20,100, 'S3'],
            ],
            vs: {
                fixed: ['0.12', '0.07', '0.17'],
                random: ['0.11', '0.03', '0.24']
            }
        }

        var vals = metajs.metaprop(d.rs, { sm: 'PFT' });
        assert.deepEqual(
            [
                vals.fixed.SM.toFixed(2), 
                vals.fixed.SM_lower.toFixed(2), 
                vals.fixed.SM_upper.toFixed(2)
            ],
            d.vs.fixed
        );
    });

    it('simple case, random incd/PFT by Inverse should be 0.11 (0.03, 0.24)', () => {
        var d = {
            rs: [
                [2, 20,  'S1'],
                [5, 90,  'S2'],
                [20,100, 'S3'],
            ],
            vs: {
                fixed: ['0.12', '0.07', '0.17'],
                random: ['0.11', '0.03', '0.24']
            }
        }

        var vals = metajs.metaprop(d.rs, { sm: 'PFT' });
        assert.deepEqual(
            [
                vals.random.SM.toFixed(2), 
                vals.random.SM_lower.toFixed(2), 
                vals.random.SM_upper.toFixed(2)
            ],
            d.vs.random
        );
    });
});


describe('testing metajs.metabin functions', () => {
    it('simple case, fixed OR by MH should be 1.697 (0.999, 2.883)', () => {
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

    it('simple case, fixed RR by MH should be 1.63 (1.00, 2.66)', () => {
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

    it('simple with 0 event, fixed OR by MH should be 1.13 (0.64, 2.00)', () => {
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

    it('simple with 0 event, fixed RR by MH should be 1.12 (0.66, 1.88)', () => {
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


    it('simple case, random OR by MH should be 2.39 (0.50, 11.45)', () => {
        var d = {
            rs: [
                [12,393,2, 396, 'S1', 'T','C'],
                [24,230,24,281, 'S2', 'T','C'],
            ],
            vs: {
                random: ['2.39', '0.50', '11.45']
            }
        };

        var vals = metajs.metabin(d.rs, { sm: 'OR' });
        assert.deepEqual(
            [
                vals.random.SM.toFixed(2), 
                vals.random.SM_lower.toFixed(2), 
                vals.random.SM_upper.toFixed(2)
            ],
            d.vs.random
        );
    });


    it('simple case, random RR by MH should be 2.34 (0.49, 11.23)', () => {
        var d = {
            rs: [
                [12,393,2, 396, 'S1', 'T','C'],
                [24,230,24,281, 'S2', 'T','C'],
            ],
            vs: {
                random: ['2.34', '0.49', '11.23']
            }
        };

        var vals = metajs.metabin(d.rs, { sm: 'RR' });
        assert.deepEqual(
            [
                vals.random.SM.toFixed(2), 
                vals.random.SM_lower.toFixed(2), 
                vals.random.SM_upper.toFixed(2)
            ],
            d.vs.random
        );
    });
});