'use strict';

import { binom } from '../src/meta.js';
import assert from 'assert'

describe('testing binom function', () => {
    it('binom.test(1, 2) should be 0.5', () => {
        assert.deepEqual(
            binom.test(1, 2),
            {
                estimate: 0.5,
                lower: 0.01258,
                upper: 0.98742
            }
        );
    });
    
    it('binom.test(2, 20) should be 0.1', () => {
        assert.deepEqual(
            binom.test(2, 20),
            {
                estimate: 0.1,
                lower: 0.01235,
                upper: 0.31698
            }
        );
    });
});
