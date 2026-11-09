//! Colored (∧^a) braid closures via cabling: each a-colored strand becomes a
//! parallel strands, a colored crossing becomes the a·b cabled crossing, and a
//! thick-edge insertion B_{w_0(S_a)} is placed at the bottom of each colored block.
//! The resulting homology is the colored homology times Π_components [a]! (as
//! q-graded spaces), up to an overall q-shift fixed by normalization.

use crate::crossing::Gen;

#[derive(Clone, Debug)]
pub struct Cabled {
    pub strands: usize,
    pub word: Vec<Gen>,
    /// start object as a word in simple reflections: product of w_0 of each block
    pub start: Vec<usize>,
    /// a-values of each colored strand (bottom positions)
    pub labels: Vec<usize>,
}

/// Longest element of S_a on positions off..off+a as a reduced word.
fn w0_word(off: usize, a: usize) -> Vec<usize> {
    let mut w = vec![];
    for i in 0..a {
        for j in (0..a - 1 - i).rev() {
            let _ = j;
        }
        let _ = i;
    }
    // standard reduced word: (s_1)(s_2 s_1)(s_3 s_2 s_1)...
    for k in 1..a {
        for j in (0..k).rev() {
            w.push(off + j);
        }
    }
    w
}

pub fn cable(labels: &[usize], word: &[Gen]) -> Cabled {
    let mut lab = labels.to_vec();
    let mut out = vec![];
    for g in word {
        let i = g.i;
        let off: usize = lab[..i].iter().sum();
        let (a, b) = (lab[i], lab[i + 1]);
        // the a-block at off..off+a crosses the b-block at off+a..off+a+b.
        // Positive σ: the left (a) block goes over. As a braid: move each strand of the
        // left block rightwards past all strands of the right block.
        // cabled word: for k = a-1 down to 0: strand at position off+k moves right b steps
        for k in (0..a).rev() {
            for j in 0..b {
                out.push(Gen { i: off + k + j, positive: g.positive });
            }
        }
        lab.swap(i, i + 1);
    }
    let mut start = vec![];
    let mut off = 0;
    for &a in labels {
        start.extend(w0_word(off, a));
        off += a;
    }
    Cabled { strands: labels.iter().sum(), word: out, start, labels: labels.to_vec() }
}
