/// Sparse parity check matrix
pub struct SparseParity {
    nrows: usize,
    ncols: usize,
    max_column_weight: usize,
    /// Each row is a row of the parity check, so M total. Each element in the
    /// row is a usize indicating a column in which a 1 resides.
    rows: Vec<Vec<usize>>,
}

fn parity_check(codeword: &[u8], parity: &SparseParity) -> Vec<u8> {
    parity.rows.iter().map(|row| {
        row.iter().copied().fold(0_u8, |acc, e| acc ^ codeword[e])
    }).collect()
}

pub fn iteration(llrs: &mut Vec<Vec<f32>>, channel_llrs: &[f32], parity: &SparseParity) {
    let mut lmn_struct = vec![];
    for _ in 0..parity.max_column_weight {
        lmn_struct.push(vec![0.0_f32; parity.ncols]);
    }
    let mut lmn_cnts = vec![0; parity.ncols];
    let mut llr_cnts = vec![0; parity.ncols];

    // parity check probs (tanh rule)
    for row in &parity.rows {
        let mut _llrs = Vec::with_capacity(row.len());
        for &c in row {
            _llrs.push(llrs[llr_cnts[c]][c]);
            llr_cnts[c] += 1;
        }
        let all = _llrs.iter().map(|l| (l/2.0).tanh()).product::<f32>();
        for &c in row {
            let lmn = all / (llrs[llr_cnts[c] - 1][c]/2.0).tanh();
            lmn_struct[lmn_cnts[c]][c] = 2.0 * lmn.atanh();
            lmn_cnts[c] += 1;
        }
    }

    println!("lmn_struct:");
    for i in 0..lmn_struct.len() {
        println!("{:.5?}", lmn_struct[i]);
    }
    println!();

    // variable check probs (repeat codes)
    let mut output = vec![0.0; parity.ncols];
    for n in 0..parity.ncols {
        let all = (0..lmn_struct.len()).map(|i| lmn_struct[i][n]).sum::<f32>();
        let mut cnt = 0;
        for row in &parity.rows {
            if row.contains(&n) {
                llrs[cnt][n] = all - lmn_struct[cnt][n] + channel_llrs[n];
                output[n] = all + channel_llrs[n];
                cnt += 1;
            }
        }
    }

    let hard_output = output.iter().map(|&o| {
        if o > 0.0 {
            0
        } else {
            1
        }
    }).collect::<Vec<u8>>();
    println!("lnm:");
    for i in 0..parity.max_column_weight {
        println!("{:.5?}", llrs[i]);
    }
    println!();
    println!("output: {:.5?}", output);
    println!();
    println!("hard output: {:?}", hard_output);
    println!();
    println!("parity check: {:?}", &parity_check(&hard_output, parity));
    println!("==================================================================================");
    println!();

}

pub fn test_sparse_decoder() {
    let h = SparseParity {
        nrows: 5,
        ncols: 10,
        max_column_weight: 3,
        rows: vec![
            vec![0, 1, 2, 5, 6, 9],
            vec![0, 2, 4, 5, 7, 8],
            vec![2, 3, 4, 6, 8, 9],
            vec![1, 3, 4, 5, 7, 9],
            vec![0, 1, 3, 6, 7, 8],
        ]
    };

    // Problem setup
    let channel_probs = vec![-0.63, -0.83, -0.73, -0.04, 0.1, 0.95, -0.76, 0.66, -0.55, 0.58];
    let channel_llrs = channel_probs.iter().map(|y| -2.0 * y).collect::<Vec<f32>>();

    let niters = 3;
    let mut llrs = vec![
        channel_llrs.clone(),
        channel_llrs.clone(),
        channel_llrs.clone()
    ];
    for _ in 0..niters {
        iteration(&mut llrs, &channel_llrs, &h);
    }
}

fn main() {
    println!("H:");
    println!("[1, 1, 1, 0, 0, 1, 1, 0, 0, 1],");
    println!("[1, 0, 1, 0, 1, 1, 0, 1, 1, 0],");
    println!("[0, 0, 1, 1, 1, 0, 1, 0, 1, 1],");
    println!("[0, 1, 0, 1, 1, 1, 0, 1, 0, 1],");
    println!("[1, 1, 0, 1, 0, 0, 1, 1, 1, 0],");
    println!();

    test_sparse_decoder();
}
