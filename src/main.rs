use generalized_djikstra::*;
use plotters::prelude::*;

const TEX_DIMS: (usize, usize) = (1024, 1024);
const OUT_FILE_NAME: &str = "plot_plain.png";

fn color_pixel(d: f32, max: f32) -> RGBColor {
    //let lambda = f32::sqrt(d / max);
    let lambda = d / max;
    let c_1: (f32, f32, f32) = (0.0, 0.0, 0.0);
    let c_2: (f32, f32, f32) = (0.0, 1.0, 0.0);
    let c_3 = (
        lambda * c_1.0 + (1.0 - lambda) * c_2.0,
        lambda * c_1.1 + (1.0 - lambda) * c_2.1,
        lambda * c_1.2 + (1.0 - lambda) * c_2.2,
    );
    // let phi = |r| f32::floor(r * 256.0) as u8;
    //RGBColor(phi(c_3.0), phi(c_3.1), phi(c_3.2))
    RGBColor(0, f32::floor(256.0 / (1.0 + 0.01 * d)) as u8, 0)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(OUT_FILE_NAME, (1024, 1024)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut vec = djikstra_eikonal::solve();
    djikstra_eikonal::djikstra(&mut vec);
    let max = vec.iter().map(|&x| f32::floor(x) as u32).max().unwrap();
    let _: Vec<Result<(), Box<dyn std::error::Error>>> = root
        .split_evenly(TEX_DIMS)
        .into_iter()
        .enumerate()
        .map(|(idx, area)| {
            let y = idx / TEX_DIMS.0;
            let x = idx % TEX_DIMS.0;
            let flip = x + TEX_DIMS.0 * (TEX_DIMS.1 - 1 - y);
            let color = color_pixel(vec[flip], max as f32);
            let _ = area.fill(&color)?;
            Ok(())
        })
        .collect();
    Ok(())
}
