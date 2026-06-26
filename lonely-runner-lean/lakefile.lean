import Lake
open Lake DSL

package «lonely-runner-lean» where
  leanOptions := #[
    ⟨`autoImplicit, false⟩
  ]

@[default_target]
lean_lib LRC where
  srcDir := "LRC"

require mathlib from git
  "https://github.com/leanprover-community/mathlib4" @ "v4.16.0"
