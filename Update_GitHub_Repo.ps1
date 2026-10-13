#requires -Version 5.1
<#
.SYNOPSIS
Update the existing 6PP ryegrass GitHub repository from the current figure scripts.
.DESCRIPTION
Reads the five scripts from ProjectRoot. Installs the merged README and runner,
then commits only the listed code/documentation paths and pushes normally.
.EXAMPLE
powershell.exe -NoProfile -ExecutionPolicy Bypass -File ".\Update_GitHub_Repo.ps1"
.EXAMPLE
.\Update_GitHub_Repo.ps1 -PreviewOnly
.EXAMPLE
.\Update_GitHub_Repo.ps1 -RepositoryPath "C:\path\to\existing\checkout" -NoPush
#>
[CmdletBinding()]
param(
    [string]$ProjectRoot = 'C:\Users\soths\OneDrive - Lincoln University\Writing\NatMicro\github',
    [string]$RepositoryPath = '',
    [string]$RepositoryUrl = 'https://github.com/Sereyboth94/6PP_Ryegrass_Metatranscriptomics.git',
    [string]$Branch = '',
    [string]$CommitMessage = 'Update manuscript figure scripts and consolidate README',
    [switch]$PreviewOnly,
    [switch]$NoPush
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

function Invoke-Git {
    param(
        [Parameter(Mandatory = $true)][string[]]$GitArgs,
        [int[]]$AcceptExitCodes = @(0),
        [switch]$Quiet
    )
    # Git writes normal progress to stderr; handle the native exit code explicitly.
    $previousErrorAction = $ErrorActionPreference
    try {
        $ErrorActionPreference = 'Continue'
        $nativeOutput = @(& git @GitArgs 2>&1)
        $nativeExit = $LASTEXITCODE
    }
    finally {
        $ErrorActionPreference = $previousErrorAction
    }
    $lines = @($nativeOutput | ForEach-Object { $_.ToString() })
    if (-not $Quiet) {
        foreach ($line in $lines) { Write-Host $line }
    }
    if ($AcceptExitCodes -notcontains $nativeExit) {
        if ($Quiet) { foreach ($line in $lines) { Write-Host $line } }
        throw "Git failed (exit $nativeExit). Stop here and resolve the reported error before rerunning."
    }
    return [pscustomobject]@{ ExitCode = $nativeExit; Output = $lines }
}

function Get-FullPath {
    param([string]$Path)
    return [System.IO.Path]::GetFullPath($Path).TrimEnd([char[]]@('\', '/'))
}

function Test-SamePath {
    param([string]$First, [string]$Second)
    return [string]::Equals((Get-FullPath $First), (Get-FullPath $Second),
        [System.StringComparison]::OrdinalIgnoreCase)
}

function Get-GitHubSlug {
    param([string]$Url)
    $clean = $Url.Trim().TrimEnd('/') -replace '\.git$', ''
    if ($clean -match '^https://github\.com/([^/]+/[^/]+)$') {
        return $Matches[1].ToLowerInvariant()
    }
    if ($clean -match '^git@github\.com:([^/]+/[^/]+)$') {
        return $Matches[1].ToLowerInvariant()
    }
    if ($clean -match '^ssh://git@github\.com/([^/]+/[^/]+)$') {
        return $Matches[1].ToLowerInvariant()
    }
    throw 'Use a GitHub HTTPS URL or the equivalent GitHub SSH URL.'
}

function Copy-SelectedFile {
    param([string]$Source, [string]$Destination)
    if (Test-SamePath $Source $Destination) { return }
    $parent = Split-Path -Path $Destination -Parent
    if (-not (Test-Path -LiteralPath $parent)) {
        New-Item -ItemType Directory -Path $parent -Force | Out-Null
    }
    Copy-Item -LiteralPath $Source -Destination $Destination -Force
}

if (-not (Get-Command git -ErrorAction SilentlyContinue)) {
    throw 'Install Git for Windows, reopen PowerShell and rerun this script.'
}
if (-not (Test-Path -LiteralPath $ProjectRoot -PathType Container)) {
    throw "Project folder not found: $ProjectRoot"
}
$ProjectRoot = Get-FullPath (Resolve-Path -LiteralPath $ProjectRoot).Path
$expectedSlug = Get-GitHubSlug $RepositoryUrl

$files = [ordered]@{
    'README.md' = (Join-Path $PSScriptRoot 'README.md')
    'Run_All_Figures.R' = (Join-Path $PSScriptRoot 'Run_All_Figures.R')
    'Update_GitHub_Repo.ps1' = $PSCommandPath
    '3.plant_biomass/script/Fig.1_Generated_Script.R' = (Join-Path $ProjectRoot '3.plant_biomass/script/Fig.1_Generated_Script.R')
    '4.plant_transcriptome/script/Figure2_Plant_Transcriptome.R' = (Join-Path $ProjectRoot '4.plant_transcriptome/script/Figure2_Plant_Transcriptome.R')
    '5.microbial_taxa/script/Figure3_Taxonomy.R' = (Join-Path $ProjectRoot '5.microbial_taxa/script/Figure3_Taxonomy.R')
    '6.microbial_function/scripts/Figure4_Function.R' = (Join-Path $ProjectRoot '6.microbial_function/scripts/Figure4_Function.R')
    '8.trait_function_linkage/script/Figure5_Trait_Function_Linkage.R' = (Join-Path $ProjectRoot '8.trait_function_linkage/script/Figure5_Trait_Function_Linkage.R')
}
foreach ($entry in $files.GetEnumerator()) {
    if (-not (Test-Path -LiteralPath $entry.Value -PathType Leaf)) {
        throw "Required source file not found: $($entry.Value)"
    }
}

if ([string]::IsNullOrWhiteSpace($RepositoryPath)) {
    $probe = Invoke-Git -GitArgs @('-C', $ProjectRoot, 'rev-parse', '--show-toplevel') -AcceptExitCodes @(0, 128) -Quiet
    if ($probe.ExitCode -eq 0 -and (Test-SamePath $ProjectRoot ($probe.Output -join ''))) {
        $RepositoryPath = $ProjectRoot
    }
    else {
        $RepositoryPath = Join-Path (Split-Path -Path $ProjectRoot -Parent) '6PP_Ryegrass_Metatranscriptomics'
    }
}
$RepositoryPath = Get-FullPath $RepositoryPath

Write-Host "Source project: $ProjectRoot"
Write-Host "Repository checkout: $RepositoryPath"
Write-Host "GitHub repository: https://github.com/$expectedSlug"
Write-Host 'Selected files:'
foreach ($path in $files.Keys) { Write-Host "  $path" }
if ($PreviewOnly) {
    Write-Host 'Preview complete. No files were copied, staged, committed or pushed.'
    return
}

if (-not (Test-Path -LiteralPath $RepositoryPath)) {
    Invoke-Git -GitArgs @('clone', '--', $RepositoryUrl, $RepositoryPath) | Out-Null
}
else {
    $probe = Invoke-Git -GitArgs @('-C', $RepositoryPath, 'rev-parse', '--show-toplevel') -AcceptExitCodes @(0, 128) -Quiet
    if ($probe.ExitCode -ne 0 -or -not (Test-SamePath $RepositoryPath ($probe.Output -join ''))) {
        $existingItems = @(Get-ChildItem -LiteralPath $RepositoryPath -Force)
        if ($existingItems.Count -eq 0) {
            Invoke-Git -GitArgs @('clone', '--', $RepositoryUrl, $RepositoryPath) | Out-Null
        }
        else {
            throw 'RepositoryPath is not the root of a Git checkout. Choose an empty/new folder with -RepositoryPath.'
        }
    }
}

$fetchRemote = Invoke-Git -GitArgs @('-C', $RepositoryPath, 'remote', 'get-url', '--all', 'origin') -AcceptExitCodes @(0, 2, 128) -Quiet
if ($fetchRemote.ExitCode -ne 0) {
    Invoke-Git -GitArgs @('-C', $RepositoryPath, 'remote', 'add', 'origin', $RepositoryUrl) | Out-Null
    $fetchRemote = Invoke-Git -GitArgs @('-C', $RepositoryPath, 'remote', 'get-url', '--all', 'origin') -Quiet
}
$pushRemote = Invoke-Git -GitArgs @('-C', $RepositoryPath, 'remote', 'get-url', '--push', '--all', 'origin') -Quiet
foreach ($url in @($fetchRemote.Output) + @($pushRemote.Output)) {
    if ((Get-GitHubSlug $url) -ne $expectedSlug) {
        throw 'An origin fetch/push URL points to a different GitHub repository. No project files were changed.'
    }
}

Invoke-Git -GitArgs @('-C', $RepositoryPath, 'fetch', 'origin') | Out-Null
if ([string]::IsNullOrWhiteSpace($Branch)) {
    $remoteHead = Invoke-Git -GitArgs @('-C', $RepositoryPath, 'symbolic-ref', '--quiet', '--short', 'refs/remotes/origin/HEAD') -AcceptExitCodes @(0, 1, 128) -Quiet
    if ($remoteHead.ExitCode -eq 0) {
        $Branch = ($remoteHead.Output -join '') -replace '^origin/', ''
    }
    else {
        $head = Invoke-Git -GitArgs @('-C', $RepositoryPath, 'ls-remote', '--symref', 'origin', 'HEAD') -Quiet
        $headLines = @($head.Output | Where-Object { $_ -match '^ref:\s+refs/heads/.+\s+HEAD$' })
        if ($headLines.Count -eq 0) {
            throw 'Could not determine the existing remote default branch. Supply -Branch with its name.'
        }
        $Branch = $headLines[0] -replace '^ref:\s+refs/heads/', '' -replace '\s+HEAD$', ''
    }
}
Invoke-Git -GitArgs @('-C', $RepositoryPath, 'check-ref-format', '--branch', $Branch) -Quiet | Out-Null
$current = Invoke-Git -GitArgs @('-C', $RepositoryPath, 'symbolic-ref', '--quiet', '--short', 'HEAD') -AcceptExitCodes @(0, 1, 128) -Quiet
if ($current.ExitCode -ne 0 -or ($current.Output -join '') -ne $Branch) {
    throw "Check out branch '$Branch', or pass -Branch for the existing branch you intend to update."
}
Write-Host "Branch: $Branch"

$remoteRef = "refs/remotes/origin/$Branch"
Invoke-Git -GitArgs @('-C', $RepositoryPath, 'rev-parse', '--verify', 'HEAD') -Quiet | Out-Null
Invoke-Git -GitArgs @('-C', $RepositoryPath, 'show-ref', '--verify', '--quiet', $remoteRef) -Quiet | Out-Null
$unmerged = Invoke-Git -GitArgs @('-C', $RepositoryPath, 'ls-files', '--unmerged') -Quiet
if ($unmerged.Output.Count -gt 0) {
    throw 'Resolve the existing merge conflicts before running the update.'
}
$counts = Invoke-Git -GitArgs @('-C', $RepositoryPath, 'rev-list', '--left-right', '--count', "HEAD...$remoteRef") -Quiet
$parts = ($counts.Output -join '').Trim() -split '\s+'
if ($parts.Count -ne 2) { throw 'Could not determine local/remote history status.' }
if ([int]$parts[0] -gt 0) {
    throw 'This checkout has pre-existing unpushed or divergent commits. Reconcile them first, or choose a fresh folder with -RepositoryPath.'
}
if ([int]$parts[1] -gt 0) {
    Invoke-Git -GitArgs @('-C', $RepositoryPath, 'merge', '--ff-only', $remoteRef) | Out-Null
}

foreach ($name in @('user.name', 'user.email')) {
    $identity = Invoke-Git -GitArgs @('-C', $RepositoryPath, 'config', '--get', $name) -AcceptExitCodes @(0, 1) -Quiet
    if ($identity.ExitCode -ne 0 -or [string]::IsNullOrWhiteSpace(($identity.Output -join ''))) {
        throw "Configure git $name with your own details, then rerun. See the README."
    }
}

foreach ($entry in $files.GetEnumerator()) {
    Copy-SelectedFile -Source $entry.Value -Destination (Join-Path $RepositoryPath $entry.Key)
}
# Install the merged documentation and runner at the source project as well.
foreach ($name in @('README.md', 'Run_All_Figures.R', 'Update_GitHub_Repo.ps1')) {
    Copy-SelectedFile -Source $files[$name] -Destination (Join-Path $ProjectRoot $name)
}

$stagePaths = @($files.Keys)
$obsoleteReadmes = @(
    'Larger_Labels_README.md',
    'Manuscript_Figures_Larger_Labels/Larger_Labels_README.md'
)
foreach ($path in $obsoleteReadmes) {
    $tracked = Invoke-Git -GitArgs @('-C', $RepositoryPath, 'ls-files', '--error-unmatch', '--', $path) -AcceptExitCodes @(0, 1) -Quiet
    if ($tracked.ExitCode -eq 0) {
        $oldFile = Join-Path $RepositoryPath $path
        if (Test-Path -LiteralPath $oldFile -PathType Leaf) {
            Remove-Item -LiteralPath $oldFile -Force
        }
        $stagePaths += $path
    }
}

Invoke-Git -GitArgs (@('-C', $RepositoryPath, 'add', '-A', '--') + $stagePaths) | Out-Null
$changed = Invoke-Git -GitArgs (@('-C', $RepositoryPath, 'diff', '--cached', '--quiet', '--') + $stagePaths) -AcceptExitCodes @(0, 1) -Quiet
if ($changed.ExitCode -eq 0) {
    Write-Host 'The selected repository files already match the current versions. No commit was needed.'
    return
}
Invoke-Git -GitArgs (@('-C', $RepositoryPath, 'diff', '--cached', '--stat', '--') + $stagePaths) | Out-Null
# Commit only these paths even if unrelated files were already staged.
Invoke-Git -GitArgs (@('-C', $RepositoryPath, 'commit', '--only', '-m', $CommitMessage, '--') + $stagePaths) | Out-Null
$sha = (Invoke-Git -GitArgs @('-C', $RepositoryPath, 'rev-parse', '--short', 'HEAD') -Quiet).Output -join ''
if ($NoPush) {
    Write-Host "Created local commit $sha. To push it:"
    Write-Host ('git -C "{0}" push --set-upstream origin "HEAD:refs/heads/{1}"' -f $RepositoryPath, $Branch)
    return
}

try {
    Invoke-Git -GitArgs @('-C', $RepositoryPath, 'push', '--set-upstream', 'origin', "HEAD:refs/heads/$Branch") | Out-Null
}
catch {
    Write-Host "Commit $sha remains in $RepositoryPath. The push did not complete; review the Git error above."
    throw
}
Write-Host "Updated https://github.com/$expectedSlug on branch $Branch (commit $sha)."

