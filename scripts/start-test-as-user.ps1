$ErrorActionPreference = "Stop"

( Enable-PSRemoting -Force -SkipNetworkProfileCheck ) | out-null
(winrm set winrm/config/service '@{AllowUnencrypted="true"}' ) | out-null
(winrm set winrm/config/service/auth '@{Basic="true"}' ) | out-null
(winrm set winrm/config/client '@{AllowUnencrypted="true"}' ) | out-null

# create credentials with no password
( $Credentials = [System.Management.Automation.PSCredential]::new("test user",[System.Security.SecureString]::new()) ) | out-null

# Set session options with NoEncryption
( $opt = New-PSSessionOption -SkipCACheck -SkipCNCheck -SkipRevocationCheck -NoEncryption ) | out-null

# Create a new session (login)

# Run the bash command as 'test user'
$env:ret = $false
Write-Output "Starting run-travis-windows.sh as user 'test user'"
Invoke-Command -Session $session -ScriptBlock {
    $script:run_script="/c/users/$env:UserName/hnn/scripts/run-travis-windows.sh"
    & "C:\Program Files\Git\bin\bash.exe" "$script:run_script"
    Write-Output "Return status: $?"
    $env:ret=$?
}

if (!$env:ret) {
    $host.SetShouldExit(-1)
    throw
}
else {
    exit 0
}
