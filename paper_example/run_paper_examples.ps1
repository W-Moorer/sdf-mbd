param(
    [string]$Configuration = "Release"
)

$ErrorActionPreference = "Stop"
$Root = Resolve-Path (Join-Path $PSScriptRoot "..")

cmake --build (Join-Path $Root "build") --config $Configuration --target demo_CH_paper_examples_openvdb
& (Join-Path $Root "build\bin\$Configuration\demo_CH_paper_examples_openvdb.exe")
python (Join-Path $PSScriptRoot "plot_paper_examples.py") --project-root $Root
python (Join-Path $PSScriptRoot "check_paper_examples.py") --project-root $Root
