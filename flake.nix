{
  description = "A Nix-flake-based R package development environment";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";

  # Track arf on main; run `nix flake update arf` (or `nix flake update`) to upgrade
  inputs.arf = {
    url = "github:eitsupi/arf";
    flake = false;
  };

  outputs =
    { self, ... }@inputs:
    let
      lib = inputs.nixpkgs.lib;
      supportedSystems = [
        "x86_64-linux"
        "aarch64-linux"
        "x86_64-darwin"
        "aarch64-darwin"
      ];
      forEachSupportedSystem =
        f:
        inputs.nixpkgs.lib.genAttrs supportedSystems (
          system:
          f {
            pkgs = import inputs.nixpkgs {
              inherit system;
              config.allowBroken = true;
              overlays = [ inputs.self.overlays.default ];
            };
          }
        );
    in
    {
      overlays.default = final: prev: rec {
        # Build arf (modern Rust-based R console) from the flake input.
        # To upgrade: run `nix flake update arf` (or `nix flake update`).
        # outputHashes only needs updating if arf changes its crossterm git pin.
        arf = final.rustPlatform.buildRustPackage {
          pname = "arf";
          version = inputs.arf.shortRev or "unstable";

          src = inputs.arf;

          cargoLock = {
            lockFile = "${inputs.arf}/Cargo.lock";
            outputHashes = {
              "crossterm-0.29.0" = "sha256-G57NGBvfZtedKQjwQMoxz1JSVH8LAPlCBeSv+DE8HiM=";
            };
          };

          # Two cd/tilde tests fail in the Nix sandbox (no $HOME), skip them
          doCheck = false;

          buildInputs = with final; lib.optionals stdenv.isDarwin [ darwin.apple_sdk.frameworks.Security ];
          nativeBuildInputs = with final; [ pkg-config ];

          meta = {
            description = "A modern Rust-based R console with fuzzy history, tree-sitter highlighting, and vi/emacs modes";
            homepage = "https://github.com/eitsupi/arf";
            license = lib.licenses.mit;
            mainProgram = "arf";
          };
        };

        # ==============================================================================
        # YOUR PACKAGE'S DEPENDENCIES (from DESCRIPTION Imports:)
        # ==============================================================================
        runtimeDeps = with final.rPackages; [
          minpack_lm
          ggplot2
          tidyr
          broom
          purrr
        ];

        # If you need packages from GitHub (not on CRAN), add them here
        # 1. Add a flake input:  inputs.myPkg = { url = "github:owner/repo"; flake = false; };
        # 2. Build it in the overlay:
        #    myPkg = final.rPackages.buildRPackage {
        #      name = "myPkg";
        #      src = inputs.myPkg;
        #      propagatedBuildInputs = with final.rPackages; [ dep1 dep2 ];
        #    };
        # 3. Reference it here:
        githubDeps = [
          # myPkg
        ];

        # ==============================================================================
        # DEVELOPMENT PACKAGES
        # ==============================================================================
        devPackages = with final.rPackages; [
          # Package development tools
          devtools
          roxygen2
          testthat
          usethis
          pkgdown
          rcmdcheck
          urlchecker

          # IDE support (R.nvim / LSP)
          httpgd
          data_table # view_df save_fun uses data.table::fwrite

          # Suggested packages (used in tests/examples)
          eemR

          # Uncomment if your package has vignettes
          # knitr
          # rmarkdown
        ];

        # Combine: runtime deps + dev tools
        rPackageList = runtimeDeps ++ githubDeps ++ devPackages;

        # ==============================================================================
        # WRAP R WITH ALL PACKAGES
        # ==============================================================================
        wrappedR = final.rWrapper.override { packages = rPackageList; };
      };

      devShells = forEachSupportedSystem (
        { pkgs }:
        {
          default = pkgs.mkShell {
            packages = with pkgs; [
              wrappedR # R with packages for LSP
              arf # modern Rust-based R console
              jarl # fast R linter (from nixpkgs)
              qpdf # PDF compression checks
              pandoc
            ];

            shellHook = ''
              export R_HOME=$(R RHOME)
              export R_LIBS_SITE=$(strings "$(command -v R)" | grep -oP '/nix/store/[^:]+/library' | sort -u | paste -sd: -)
              export R_LIBS_USER="$PWD/.r-libs"
              mkdir -p "$R_LIBS_USER"
              export R_QPDF="${pkgs.qpdf}/bin/qpdf"
            '';
          };
        }
      );
    };
}
