//! Build and persist GrTable(r, b) tables (M4c B2): `grwarm r:b [r:b ...]`.
//! Each table is loaded from the on-disk cache if present and valid, else computed and written
//! (see `symfun::GrTable::load_or_compute`). Prints the wall time and the cache file.
use krhom::algebra::symfun::GrTable;
fn main() {
    let dir = GrTable::cache_dir().expect("KRHOM_GRCACHE=off: nothing to do");
    for a in std::env::args().skip(1) {
        let v: Vec<usize> = a.split(':').map(|x| x.parse().expect("r:b")).collect();
        let (r, b) = (v[0], v[1]);
        let t0 = std::time::Instant::now();
        let hit = GrTable::load(&GrTable::cache_path(&dir, r, b), r, b).is_some();
        let t = GrTable::get(r, b);
        if !GrTable::cache_path(&dir, r, b).exists() {
            t.save(&dir).expect("save");
        }
        println!("GrTable({}, {}) fsize {} {} in {:.1} s -> {}", r, b, t.fsize(), if hit { "loaded" } else { "computed" }, t0.elapsed().as_secs_f64(), GrTable::cache_path(&dir, r, b).display());
    }
}
