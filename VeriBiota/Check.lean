import Lean

open IO System

private def runLean (path : FilePath) : IO UInt32 := do
  let child ← Process.spawn {cmd := "lean", args := # ["--check", path.toString], stdout := .inherit, stderr := .piped}
  let exit ← child.wait
  if exit ≠ 0 then
    let err ← child.stderr.readToEnd
    if err.trim.isEmpty then
      IO.eprintln s!"lean --check {path} failed with exit code {exit}."
    else
      IO.eprintln err
  pure exit

@[export main]
def main (args : List String) : IO UInt32 := do
  match args with
  | [] =>
      IO.eprintln "usage: lake exe veribiota-check <Lean file>"
      pure 1
  | file :: _ => do
      let path := FilePath.mk file
      let exists ← path.pathExists
      if !exists then
        IO.eprintln s!"Lean file '{path}' not found."
        pure 1
      else
        runLean path
