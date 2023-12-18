
use std::collections::HashMap;

use petgraph::{prelude::*, algo::toposort};
use crate::alignment::{*, self};

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

        self.add_more_read();
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

        let align_matrix_forward = self.alignment(read);
        let alin_matrix_reverse_complement = 0;

    }

    fn alignment(&self, read: &[u8]) -> AlignMatrix{

        let mut nodes = toposort(&self.graph, None).unwrap();
        assert!(nodes[0] == self.begin_node_idx);
        assert!(*nodes.last().unwrap() == self.end_node_idx);
        
        // remove the begin and end node
        nodes.pop();
        nodes.remove(0);

        let node_idx_to_col_idx = nodes
            .iter()
            .enumerate()
            .map(|v| (*v.1, v.0))
            .collect::<HashMap<_, _>>();

        let matrix_rows = read.len() + 1;
        let matrix_cols = nodes.len() + 1;
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
                let node_idx = nodes[col - 1];
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