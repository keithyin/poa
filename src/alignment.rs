pub struct AlignConfig {
    pub align_mode: AlignMode,
    pub align_params: AlignParams
}

impl AlignConfig {
    pub fn new(align_mode: AlignMode, align_params: AlignParams) -> Self {
        AlignConfig { align_mode, align_params}
    }

    pub fn default() -> Self {
        AlignConfig::new(AlignMode::LOCAL, AlignParams::default())
    }
}

pub enum AlignMode {
    GLOBAL,
    SEMIGLOBAL,
    LOCAL
}


pub struct AlignParams {
    pub match_score: i32,
    pub mismatch_score: i32,
    pub insertion_score: i32,
    pub deletion_score: i32,
}

impl AlignParams {
    pub fn new(match_score: i32, mismatch_score: i32, insertion_score: i32, deletion_score: i32) -> Self {
        AlignParams{match_score, mismatch_score, insertion_score, deletion_score}
    }

    pub fn default() -> Self {
        AlignParams::new(3, -5, -4, -4)
    }
}

#[derive(Clone, Copy)]
pub enum TransMode {
    START,
    MATCH,
    INSERT,
    DELETE,
    END
}

#[derive(Clone, Copy)]
pub struct AlignPosition {
    pub score: i32,
    pub pre_trans_mode: TransMode, // how to trans to this position!
    pub pre_column: i32
}

impl AlignPosition {
    pub fn new(score: i32, pre_trans_mode: TransMode, pre_colunm: i32) -> Self {
        AlignPosition{
            score, 
            pre_trans_mode, 
            pre_column: pre_colunm
        }
    }

    pub fn default() -> Self {
        AlignPosition::new(0, TransMode::START, -1)
    }
}

pub struct AlignMatrix {
    matrix: Vec<Vec<Option<AlignPosition>>>,
    max_score: i32,
    max_score_positions: Vec<(usize, usize)>
}

impl AlignMatrix {
    
    pub fn new(num_rows: usize, num_cols: usize) -> Self {
        AlignMatrix{
            matrix: vec![vec![None; num_cols]; num_rows],
            max_score: i32::MIN,
            max_score_positions: vec![],
        }
    }

    pub fn set(&mut self, row: usize, col: usize, align_pos: AlignPosition) {
        if align_pos.score > self.max_score {
            self.max_score = align_pos.score;
            self.max_score_positions.clear();
        }
        
        if align_pos.score == self.max_score {
            self.max_score_positions.push((row, col));
        }

        self.matrix[row][col] = Some(align_pos);
    }

    pub fn get(&self, row: usize, col: usize) -> Option<AlignPosition> {
        return self.matrix[row][col];
    }

}


///    '' A T C G
/// '' 0
/// A  0
/// T
/// C
/// G
/// global matrix[i, j] = max(
///                            score(s1[i-1], s2[j-1]) + matrix[i-1, j-1],
///                            insertion + matrix[i-1, j],       to bottom                     
///                            deletion + matrix[i, j-1]         to right
///  )
/// local matrix[i, j] = max(   0, 
///                            score(s1[i-1], s2[j-1]) + matrix[i-1, j-1],
///                            insertion + matrix[i-1, j],       to bottom                     
///                            deletion + matrix[i, j-1]         to right
///  )
/// s1 is query, s2 is target.
pub fn align(s1: &[u8], s2: &[u8], align_cfg: &AlignConfig) -> AlignMatrix{
    let matrix_rows = s1.len() + 1;
    let matrix_cols = s2.len() + 1;
    let mut align_matrix = AlignMatrix::new(matrix_rows, matrix_cols);
    
    for row in 0..matrix_rows {
        match align_cfg.align_mode {
            AlignMode::GLOBAL => align_matrix.set(row, 0, 
                AlignPosition::new(align_cfg.align_params.insertion_score * (row as i32), TransMode::INSERT, -1)),
            AlignMode::LOCAL => align_matrix.set(row, 0, AlignPosition::new(0, TransMode::START, -1)),
            _ => panic!("invalid align_mode"),
        }
    }

    for col in 0..matrix_cols {
        match align_cfg.align_mode {
            AlignMode::GLOBAL => align_matrix.set(0, col, 
                AlignPosition::new(align_cfg.align_params.insertion_score * (col as i32), TransMode::DELETE, -1)),
            AlignMode::LOCAL => align_matrix.set(0, col, AlignPosition::new(0, TransMode::START, -1)),
            _ => panic!("invalid align_mode"),
        }
    }

    for col in 1..matrix_cols {
        for row in 1..matrix_rows {
            
            let mut align_position =  match align_cfg.align_mode{
                AlignMode::GLOBAL => AlignPosition::new(i32::MIN, TransMode::START, -1),
                AlignMode::LOCAL => AlignPosition::new(0, TransMode::START, -1),
                _ => panic!("invalid align_mode"),
            };

            // match
            let mut score = if s1[row-1] == s2[col-1] {align_cfg.align_params.match_score} else {align_cfg.align_params.mismatch_score} 
                + align_matrix.get(row -1, col-1).unwrap().score;
            
            if score > align_position.score {
                align_position.score = score;
                align_position.pre_trans_mode = TransMode::MATCH;
            }
            
            // deletion
            score = align_cfg.align_params.deletion_score + align_matrix.get(row, col-1).unwrap().score;
            if score > align_position.score {
                align_position.score = score;
                align_position.pre_trans_mode = TransMode::DELETE;
            }

            // insertion
            score = align_cfg.align_params.insertion_score + align_matrix.get(row-1, col).unwrap().score;
            if score > align_position.score {
                align_position.score = score;
                align_position.pre_trans_mode = TransMode::INSERT;
            }

            align_matrix.set(row, col, align_position);
        }
    }
    align_matrix
}
