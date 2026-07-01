# Parse THERM.DAT and generate mechanism_data.ipp
param(
    [string]$ThermPath = "",
    [string]$OutPath = ""
)

if ($ThermPath -eq "") {
    $ThermPath = Join-Path $PSScriptRoot "..\..\detailed-CH4-Air\THERM.DAT"
}
if ($OutPath -eq "") {
    $OutPath = Join-Path $PSScriptRoot "mechanism_data.ipp"
}

Write-Host "ThermPath: $ThermPath"
Write-Host "OutPath: $OutPath"

if (-not (Test-Path $ThermPath)) {
    Write-Error "THERM.DAT not found at: $ThermPath"
    exit 1
}

$therm = Get-Content $ThermPath -Encoding UTF8 -Raw

# Species we need
$targets = @(
    'H', 'H2', 'O', 'O2', 'OH', 'HO2', 'H2O', 'H2O2',
    'CO', 'CO2', 'HCO', 'CH3', 'CH4', 'CH2O', 'CH2', 'CH3O',
    'CH2OH', 'CH', 'C', 'SCH2', 'CH3O2', 'CH3O2H',
    'C2H6', 'C2H5', 'C2H4', 'C2H3', 'C2H2', 'C2H', 'CH2CO', 'HCCO',
    'CH3OH', 'CH3HCO', 'CH3CO', 'CH2HCO', 'C2H4O', 'C2H5O', 'C2H5O2',
    'C2H5O2H', 'C2', 'C2O', 'CH3CO3', 'CH3CO3H',
    'N2'
)

# Parse species blocks
$lines = $therm -split '\r?\n'
$specData = @{}
$i = 0
while ($i -lt $lines.Count) {
    $line = $lines[$i]
    if ($line.Length -ge 45 -and $line -match '^\s*(\S+)') {
        $name = $Matches[1]
        if ($name -eq 'THERMO' -or $name -eq 'END') { $i++; continue }
        if ($targets -contains $name -and -not $specData.ContainsKey($name)) {
            $block = @($line)
            for ($k = 1; $k -lt 4 -and ($i + $k) -lt $lines.Count; $k++) {
                $block += $lines[$i + $k]
            }
            $specData[$name] = $block
        }
    }
    $i++
}

Write-Host "Found $($specData.Count)/$($targets.Count) species"

# Molecular weights from chem.out
$mw = @{
    'H'='1.00797e-3'; 'H2'='2.01594e-3'; 'O'='15.99940e-3'; 'O2'='31.99880e-3';
    'OH'='17.00737e-3'; 'HO2'='33.00677e-3'; 'H2O'='18.01534e-3'; 'H2O2'='34.01474e-3';
    'CO'='28.01055e-3'; 'CO2'='44.00995e-3'; 'HCO'='29.01852e-3'; 'CH3'='15.03506e-3';
    'CH4'='16.04303e-3'; 'CH2O'='30.02649e-3'; 'CH2'='14.02709e-3'; 'CH3O'='31.03446e-3';
    'CH2OH'='31.03446e-3'; 'CH'='13.01912e-3'; 'C'='12.01115e-3'; 'SCH2'='14.02709e-3';
    'CH3O2'='47.03386e-3'; 'CH3O2H'='48.04183e-3';
    'C2H6'='30.07012e-3'; 'C2H5'='29.06215e-3'; 'C2H4'='28.05418e-3'; 'C2H3'='27.04621e-3';
    'C2H2'='26.03824e-3'; 'C2H'='25.03027e-3'; 'CH2CO'='42.03764e-3'; 'HCCO'='41.02967e-3';
    'CH3OH'='32.04243e-3'; 'CH3HCO'='44.05358e-3'; 'CH3CO'='43.04561e-3'; 'CH2HCO'='43.04561e-3';
    'C2H4O'='44.05358e-3'; 'C2H5O'='45.06155e-3'; 'C2H5O2'='61.06095e-3'; 'C2H5O2H'='62.06892e-3';
    'C2'='24.02230e-3'; 'C2O'='40.02170e-3'; 'CH3CO3'='75.04441e-3'; 'CH3CO3H'='76.05238e-3';
    'N2'='28.01340e-3'
}

$sb = [System.Text.StringBuilder]::new()

[void]$sb.AppendLine("/* Auto-generated from THERM.DAT */")
[void]$sb.AppendLine("// Species: $($targets.Count) C/H/O species + N2")
[void]$sb.AppendLine("")

[void]$sb.AppendLine("const char* SPECIES_NAMES_HOST[] = {")
foreach ($name in $targets) { [void]$sb.AppendLine("    `"$name`",") }
[void]$sb.AppendLine("};")
[void]$sb.AppendLine("")

[void]$sb.AppendLine("real MOL_WT_HOST[] = {")
foreach ($name in $targets) { [void]$sb.AppendLine("    $($mw[$name])f,") }
[void]$sb.AppendLine("};")
[void]$sb.AppendLine("")

[void]$sb.AppendLine("real NASA_TMID_HOST[] = {")
foreach ($name in $targets) {
    if ($specData.ContainsKey($name)) {
        $line1 = $specData[$name][0]
        if ($line1.Length -ge 65) {
            $tStr = $line1.Substring(45, [Math]::Min(30, $line1.Length - 45))
            $tParts = $tStr -split '\s+' | Where-Object { $_ }
            if ($tParts.Count -ge 3) {
                [void]$sb.AppendLine("    $($tParts[2])f,  // $name")
            } else {
                [void]$sb.AppendLine("    1000.0f,  // $name (default)")
            }
        } else {
            [void]$sb.AppendLine("    1000.0f,  // $name (default)")
        }
    } else {
        [void]$sb.AppendLine("    1000.0f,  // $name (missing)")
    }
}
[void]$sb.AppendLine("};")
[void]$sb.AppendLine("")

[void]$sb.AppendLine("real NASA_COEFFS_HOST[][14] = {")
foreach ($name in $targets) {
    if ($specData.ContainsKey($name)) {
        $block = $specData[$name]
        $coeffText = $block[1] + $block[2] + $block[3]
        $coeffText = $coeffText -replace '\s+', ' '
        $coeffMatches = [regex]::Matches($coeffText, '([+-]?\d+\.\d+E[\s+\-]?\d+)')
        $coeffs = @()
        foreach ($m in $coeffMatches) { 
            $val = $m.Value -replace '\s', ''
            $val = $val -replace 'E(?![\+\-])', 'E+'
            $coeffs += [double]::Parse($val) 
        }
        
        if ($coeffs.Count -ge 14) {
            $c14 = $coeffs[0..13]
            $cStr = ($c14 | ForEach-Object { "$($_.ToString('E8'))f" }) -join ", "
            [void]$sb.AppendLine("    { $cStr },  // $name")
        } elseif ($coeffs.Count -ge 7) {
            $c7 = $coeffs[0..6]
            $cStr = ($c7 | ForEach-Object { "$($_.ToString('E8'))f" }) -join ", "
            [void]$sb.AppendLine("    { $cStr, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f },  // $name (single range)")
        } else {
            [void]$sb.AppendLine("    { 0.0f },  // $name (ERROR)")
        }
    } else {
        [void]$sb.AppendLine("    { 0.0f },  // $name (missing)")
    }
}
[void]$sb.AppendLine("};")
[void]$sb.AppendLine("")

[void]$sb.AppendLine("void copyMechanismToDevice() {")
[void]$sb.AppendLine("    cudaMemcpyToSymbol(MOL_WT, MOL_WT_HOST, sizeof(real) * NUM_SPECIES);")
[void]$sb.AppendLine("    cudaMemcpyToSymbol(NASA_TMID, NASA_TMID_HOST, sizeof(real) * NUM_SPECIES);")
[void]$sb.AppendLine("    cudaMemcpyToSymbol(NASA_COEFFS, NASA_COEFFS_HOST, sizeof(real) * NUM_SPECIES * 14);")
[void]$sb.AppendLine("}")

$sb.ToString() | Out-File -FilePath $OutPath -Encoding UTF8
Write-Host "Generated $OutPath"