
use std::collections::HashMap;

use petgraph::{prelude::*, algo::toposort};
use crate::{alignment::{*, self}, dna_utils::reverse_complement};

pub struct PartialOrderAlignment {

    graph: Graph<PoaNode, PoaEdge>,
    begin_node_idx: NodeIndex,
    end_node_idx: NodeIndex,
    num_reads: usize,

    align_cfg: AlignConfig,
    score_threshold: f32,
}


impl PartialOrderAlignment {
    
    pub fn new(align_config: AlignConfig, score_threshold: f32) -> Self{
        let mut graph = Graph::new();
        let begin_node_idx = graph.add_node(PoaNode::new('-' as u8, 0));
        let end_node_idx = graph.add_node(PoaNode::new('-' as u8, 0));
        PartialOrderAlignment{
            graph, 
            begin_node_idx, 
            end_node_idx, 
            num_reads: 0, 
            align_cfg: align_config,
            score_threshold: score_threshold}
    }

    pub fn add_read(&mut self, read: &[u8]) {
        if self.num_reads == 0 {
            self.add_first_read(read);
            return;
        }

        self.add_more_read(read);
    }

    fn add_first_read(&mut self, read: &[u8]) {
        assert!(read.len() > 0 && self.num_reads == 0, "read.len={}, self.num_reads={}", read.len(), self.num_reads);
        let mut pre_node_idx = self.begin_node_idx;
        for base in read {
            let cur_node_idx = self.graph.add_node(PoaNode::new(*base, 1));
            self.graph.add_edge(pre_node_idx, cur_node_idx, PoaEdge::new(0));
            pre_node_idx = cur_node_idx;
        }
        self.graph.add_edge(pre_node_idx, self.end_node_idx, PoaEdge::new(0));
        self.num_reads += 1;
    }

    fn add_more_read(&mut self, read: &[u8]) {
        let mut target = toposort(&self.graph, None).unwrap();
        assert!(target[0] == self.begin_node_idx);
        assert!(*target.last().unwrap() == self.end_node_idx);
        target.pop();
        target.remove(0);

        let align_matrix_forward = self.alignment(&target, read);
        let reverse_seq = reverse_complement(read);
        let align_matrix_reverse_complement = self.alignment(&target, &reverse_seq);
        if align_matrix_forward.get_max_score() > align_matrix_reverse_complement.get_max_score() {
            self.commit_add(&target, read, &align_matrix_forward);
        } else {
            self.commit_add(&target, &reverse_seq, &align_matrix_reverse_complement);
        }
        self.num_reads += 1;

    }

    fn alignment(&self, target: &Vec<NodeIndex>, read: &[u8]) -> AlignMatrix{
        let node_idx_to_col_idx = target
            .iter()
            .enumerate()
            .map(|v| (*v.1, v.0))
            .collect::<HashMap<_, _>>();

        let matrix_rows = read.len() + 1;
        let matrix_cols = target.len() + 1;
        let mut align_matrix = AlignMatrix::new(matrix_rows, matrix_cols);

        for row in 0..matrix_rows {
            match self.align_cfg.align_mode {
                AlignMode::GLOBAL => align_matrix.set(row, 0, 
                    AlignPosition::new(self.align_cfg.align_params.insertion_score * (row as i32), TransMode::INSERT, -1)),
                AlignMode::LOCAL => align_matrix.set(row, 0, AlignPosition::new(0, TransMode::START, -1)),
                _ => panic!("invalid align_mode"),
            }
        }
    
        for col in 0..matrix_cols {
            match self.align_cfg.align_mode {
                AlignMode::GLOBAL => align_matrix.set(0, col, 
                    AlignPosition::new(self.align_cfg.align_params.insertion_score * (col as i32), TransMode::DELETE, -1)),
                AlignMode::LOCAL => align_matrix.set(0, col, AlignPosition::new(0, TransMode::START, -1)),
                _ => panic!("invalid align_mode"),
            }
        }
        
        for row in 1..matrix_rows {
            for col in 1..matrix_cols {
                let node_idx = target[col - 1];
                let graph_node = self.graph.node_weight(node_idx).unwrap();

                let mut align_position = match self.align_cfg.align_mode {
                    AlignMode::GLOBAL => AlignPosition::new(i32::MIN, TransMode::START, -1),
                    AlignMode::LOCAL => AlignPosition::new(0, TransMode::START, -1),
                    _ => panic!("invalid AlignMode"),
                };

                for edge in self.graph.edges_directed(node_idx, Incoming) {
                    let pre_node_idx = edge.source();
                    let pre_col = node_idx_to_col_idx[&pre_node_idx] + 1;

                    // match
                    let mut score = if graph_node.base == read[row-1] {
                        self.align_cfg.align_params.match_score
                    } else {
                        self.align_cfg.align_params.match_score
                    } + align_matrix.get(row-1, pre_col).unwrap().score;

                    if score > align_position.score {
                        align_position.score = score;
                        align_position.pre_column = pre_col as i32;
                        align_position.pre_trans_mode = TransMode::MATCH;
                    }
                    
                    // deletion
                    score = self.align_cfg.align_params.insertion_score + align_matrix.get(row, pre_col).unwrap().score;
                    if score > align_position.score {
                        align_position.score = score;
                        align_position.pre_column = pre_col as i32;
                        align_position.pre_trans_mode = TransMode::DELETE;
                    }

                    // insertion
                    score = self.align_cfg.align_params.insertion_score + align_matrix.get(row - 1, col).unwrap().score;
                    if score > align_position.score {
                        align_position.score = score;
                        align_position.pre_column = col as i32;
                        align_position.pre_trans_mode = TransMode::INSERT;
                    }
                }
                align_matrix.set(row, col, align_position);
            }
        }
        align_matrix
    }

    /// 对齐的部分接上
    /// 没有对齐的部分，接到头尾上
    /// align_matrix: row is read, target is col
    fn commit_add(&mut self, target: &Vec<NodeIndex>, read: &[u8], align_matrix: &AlignMatrix) {
        let max_score_position = align_matrix.get_max_score_positions().last().unwrap();
        
        let read_end = max_score_position.0 - 1;
        let target_end = max_score_position.1 - 1;

        let mut read_cursor = read.len() - 1;
        let mut target_cursor = target_end;
        // read_cursor 对应的是 read 的索引，对于 matrix 的索引 需要+1
        let mut next_node = self.end_node_idx;
        
        while read_cursor > read_end {
            let new_node = self.graph.add_node(PoaNode::new(read[read_cursor], 1));
            self.graph.add_edge(new_node, next_node, PoaEdge::new(0));
            next_node = new_node;
            read_cursor -= 1;
        }

        loop {
            let align_pos = align_matrix.get(read_cursor+1, target_cursor+1).unwrap(); 

            let skip = match align_pos.pre_trans_mode {
                TransMode::START => {
                    true
                },
                TransMode::MATCH => {
                    let mut target_node = self.graph.node_weight_mut(target[target_cursor]).unwrap();
                    if target_node.base == read[read_cursor] {
                        target_node.num_reads += 1;
                        next_node = target[target_cursor];
                    } else {
                        let new_node = self.graph.add_node(PoaNode::new(read[read_cursor], 1));
                        self.graph.add_edge(new_node, next_node, PoaEdge::new(0));
                        next_node = new_node;
                    }
                    read_cursor -= 1;
                    target_cursor -= 1;
                    false
                },
                TransMode::DELETE => {
                    target_cursor -= 1;
                    false
                },
                TransMode::INSERT => {
                    let new_node = self.graph.add_node(PoaNode::new(read[read_cursor], 1));
                    self.graph.add_edge(new_node, next_node, PoaEdge::new(0));
                    next_node = new_node;
                    read_cursor -= 1;
                    false
                },
                _ => panic!("invalid"),
            };

            if skip {
                break;
            }
        }

        loop {
            let new_node = self.graph.add_node(PoaNode::new(read[read_cursor], 1));
            self.graph.add_edge(new_node, next_node, PoaEdge::new(0));
            next_node = new_node;
            if read_cursor == 0 {
                break;
            }
            read_cursor -= 1;
        }

        self.graph.add_edge(self.begin_node_idx, next_node, PoaEdge::new(0));

        

    }
}


#[derive(Debug, Clone, Copy)]
pub struct PoaNode {
    base: u8,
    num_reads: u32,
    num_spinning_reads: u32

}

impl PoaNode {
    pub fn new(base: u8, num_reads: u32) -> Self {
        PoaNode{base, num_reads, num_spinning_reads: 0}
    }
}


#[derive(Debug, Clone, Copy)]
pub struct PoaEdge {
    weight: u32
}

impl PoaEdge {
    pub fn new(weight: u32) -> Self {
        PoaEdge { weight }
    }
}



#[cfg(test)]
mod test {

    #[test]
    fn test_poa() {
        println!("hello");
    }
}