use priority_queue::PriorityQueue;
use std::cmp::Ordering;
use std::collections::HashSet;

const TEX_DIMS: (u32, u32) = (1024, 1024);
const WINDOW: ((f32, f32), (f32, f32)) = ((-3.0, 3.0), (-3.0, 3.0));
const EPSILON: f32 = 0.01;

fn idx_to_point(idx: usize) -> (i32, i32) {
    let cols = TEX_DIMS.0 as i32;
    (
        i32::rem_euclid(idx as i32, cols),
        i32::div_euclid(idx as i32, cols),
    )
}

fn point_to_idx(p: (i32, i32)) -> Option<usize> {
    if p.0 < TEX_DIMS.0 as i32 && p.0 >= 0 && p.1 < TEX_DIMS.1 as i32 && p.1 >= 0 {
        Some((p.0 + (p.1) * (TEX_DIMS.0 as i32)) as usize)
    } else {
        None
    }
}

pub fn solve() -> Vec<f32> {
    let x_init = WINDOW.0 .0 as f32;
    let x_step = (WINDOW.0 .1 - WINDOW.0 .0) / (TEX_DIMS.0 as f32);
    let y_init = WINDOW.1 .0 as f32;
    let y_step = (WINDOW.1 .1 - WINDOW.1 .0) / (TEX_DIMS.1 as f32);
    (0..(TEX_DIMS.0 * TEX_DIMS.1))
        .into_iter()
        .map(|idx| {
            let (x, y) = idx_to_point(idx as usize);
            let x_real = x_init + (x as f32) * x_step;
            let y_real = y_init + (y as f32) * y_step;
            if f32::abs(to_plot(x_real, y_real)) < EPSILON {
                0.0
            } else {
                1.0 / 0.0
            }
        })
        .collect()
}

fn to_plot(x: f32, y: f32) -> f32 {
    return x * x + y * y - (0.5 * (2.0 + f32::sin(5.0 * f32::atan2(y, x))));
}

fn val_at_point(p: (i32, i32), vals: &Vec<f32>) -> f32 {
    match point_to_idx(p) {
        Some(idx) => vals[idx],
        None => 1.0 / 0.0,
    }
}

fn eikonal_update(p: (i32, i32), vals: &Vec<f32>) -> f32 {
    let up_p = (p.0 + 1, p.1);
    let down_p = (p.0 - 1, p.1);
    let left_p = (p.0, p.1 - 1);
    let right_p = (p.0, p.1 + 1);
    let u_x = f32::min(val_at_point(up_p, vals), val_at_point(down_p, vals));
    let u_y = f32::min(val_at_point(right_p, vals), val_at_point(left_p, vals));
    let a = f32::max(u_x, u_y);
    let b = f32::min(u_x, u_y);
    if b == 1.0 / 0.0 {
        1.0 / 0.0
    } else if a - b > 1.0 {
        b + 1.0
    } else {
        0.5 * (u_x + u_y + f32::sqrt(2.0 - (a - b) * (a - b)))
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
struct CoOrd(u32, u32);

#[derive(Clone, Copy, PartialEq, Debug)]
struct MinFloat(f32);

impl Eq for MinFloat {}

impl PartialOrd for MinFloat {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        return other.0.partial_cmp(&self.0);
    }
}

impl Ord for MinFloat {
    fn cmp(&self, other: &Self) -> Ordering {
        return self.partial_cmp(other).unwrap();
    }
}

pub fn djikstra(vals: &mut Vec<f32>) {
    let mut solved_points: HashSet<(i32, i32)> = HashSet::new();
    let mut pq: PriorityQueue<(i32, i32), MinFloat> = PriorityQueue::new();
    vals.iter().enumerate().for_each(|(idx, &val)| {
        pq.push(idx_to_point(idx), MinFloat(val));
    });
    while let Some((p, _)) = pq.pop() {
        solved_points.insert(p);
        let unsolved_neibs: Vec<(i32, i32)> = [(1, 0), (-1, 0), (0, 1), (0, -1)]
            .into_iter()
            .map(|(u, v)| (p.0 + u, p.1 + v))
            .filter(|&p0| point_to_idx(p0).is_some() && !solved_points.contains(&p0))
            .collect();
        unsolved_neibs.into_iter().for_each(|p0| {
            let r0 = val_at_point(p0, vals);
            let r1 = eikonal_update(p0, vals);
            let r2 = f32::min(r0, r1);
            vals[point_to_idx(p0).unwrap()] = r2;
            let _ = pq.push(p0, MinFloat(r2));
        });
    }
}
