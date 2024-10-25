const H: [[u8; 10]; 5] = [
    [1, 1, 1, 0, 0, 1, 1, 0, 0, 1],
    [1, 0, 1, 0, 1, 1, 0, 1, 1, 0],
    [0, 0, 1, 1, 1, 0, 1, 0, 1, 1],
    [0, 1, 0, 1, 1, 1, 0, 1, 0, 1],
    [1, 1, 0, 1, 0, 0, 1, 1, 1, 0],
];


fn check_node_to_variables(llrs: &[f64]) -> f64 {
    2.0 * (llrs.iter().map(|l| (l/2.0).tanh()).product::<f64>()).atanh()
}

fn variable_to_check_nodes(llrs: &[f64], ln: f64) -> f64 {
    ln + llrs.iter().sum::<f64>()
}

fn iteration(llrs: &mut [[f64; 10]; 3], channel_llrs: &[f64]) {
    let mut lmn_struct: [[f64; 10]; 3] = [
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ];
    let mut lmn_cnts: [usize; 10] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    let mut llr_cnts: [usize; 10] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

    // parity check probs (tanh rule)
    for m in 0..5 {
        let mut _llrs = vec![];
        for i in 0..10 {
            if H[m][i] == 1 {
                _llrs.push(llrs[llr_cnts[i]][i]);
                llr_cnts[i] += 1;
            }
        }
        let all = _llrs.iter().map(|l| (l/2.0).tanh()).product::<f64>();
        for n in 0..10 {
            if H[m][n] == 1 {
                let lmn = all / (llrs[llr_cnts[n] - 1][n]/2.0).tanh();
                lmn_struct[lmn_cnts[n]][n] = 2.0 * lmn.atanh();
                lmn_cnts[n] += 1;
            }
        }
    }

    println!("lmn_struct:");
    for i in 0..3 {
        println!("{:.5?}", lmn_struct[i]);
    }
    println!();

    // variable check probs (repeat codes)
    let mut output: [f64; 10] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    for n in 0..10 {
        let all = (0..3).map(|i| lmn_struct[i][n]).sum::<f64>();
        let mut cnt = 0;
        for m in 0..5 {
            if H[m][n] == 1 {
                llrs[cnt][n] = all - lmn_struct[cnt][n] + channel_llrs[n];
                output[n] = all + channel_llrs[n];
                cnt += 1;
            }
        }
    }
    
    let hard_output: [u8; 10] = output.iter().map(|&o| {
        if o > 0.0 {
            0
        } else {
            1
        }
    }).collect::<Vec<_>>().try_into().unwrap();
    println!("lnm:");
    for i in 0..3 {
        println!("{:.5?}", llrs[i]);
    }
    println!();
    println!("output: {:.5?}", output);
    println!();
    println!("hard output: {:?}", hard_output);
    println!();
    println!("parity check: {:?}", &parity_check(hard_output));
    println!("==================================================================================");
    println!();
}

fn parity_check(codeword: [u8; 10]) -> [u8; 5] {
    let mut output: [u8; 5] = [0, 0, 0, 0, 0];
    for m in 0..5 {
        for n in 0..10 {
            if H[m][n] == 1 {
                output[m] ^= codeword[n];
            }
        }
    }
    output
}


fn main() {
    println!("H:");
    for row in H {
        println!("{:?}", row);
    }
    println!();

    // Problem setup
    let channel_probs: [f64; 10] = [-0.63, -0.83, -0.73, -0.04, 0.1, 0.95, -0.76, 0.66, -0.55, 0.58];
    let channel_llrs: [f64; 10] = channel_probs.iter().map(|y| -2.0 * y).collect::<Vec<_>>().try_into().unwrap();
    //let channel_llrs: [f64; 10] = [1.3, 1.7, 1.5, 0.08, -0.2, -1.9, 1.5, -1.3, 1.1, -1.2];

    let niters = 3;
    let mut llrs = [
        channel_llrs.clone(),
        channel_llrs.clone(),
        channel_llrs.clone()
    ];
    for _ in 0..niters {
        iteration(&mut llrs, &channel_llrs);
    }
}
