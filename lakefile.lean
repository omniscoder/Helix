import Lake
open Lake DSL

package helixVeriBiota

@[default_target]
lean_lib HelixVeriBiota where
  globs := #[`VeriBiota]

lean_exe veribiota-check where
  root := `VeriBiota.Check
