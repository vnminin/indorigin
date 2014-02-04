test_that("q probs sum approximately to transition probs", 
        {
            l0 = 1.0
            l1 = 2.3
            t = 4.76
            trans = two.state.trans.prob(l0, l1, t)
            n_max = 200
            trans_q = pcltTest(l0, l1, t, n_max)
            
            q01 = trans[1, 2]
            q01_test = sum(exp(trans_q[seq(2, n_max, by = 2)]))
            expect_equal(q01_test, q01, tolerance = 1e-6,
                    scale = q01)
            q00 = trans[1, 1]
            q00_test = sum(exp(trans_q[seq(1, n_max, by = 2)]))
            expect_equal(q00_test, q00, tolerance = 1e-6,
                    scale = q00)
        })