Enable-WindowsOptionalFeature -Online -FeatureName Microsoft-Windows-Subsystem-Linux
Invoke-WebRequest -Uri https://aka.ms/wsl-ubuntu-1804 -OutFile Ubuntu.appx -UseBasicParsing
Copy-Item .\Ubuntu.appx .\Ubuntu.zip
Expand-Archive .\Ubuntu.zip .\Ubuntu
$userenv = [System.Environment]::GetEnvironmentVariable("Path", "User"); [System.Environment]::SetEnvironmentVariable("PATH", $userenv + ";C:\Users\Travis\Ubuntu", "User")

$Credentials = [System.Management.Automation.PSCredential]::new("travis",[System.Security.SecureString]::new())
$opt = New-PSSessionOption -SkipCACheck -SkipCNCheck -SkipRevocationCheck -NoEncryption
$session = New-PSSession -ComputerName localhost -SessionOption $opt -Authentication Basic -Credential $Credentials
Invoke-Command -Session $session -ScriptBlock {Add-AppxPackage .\Ubuntu.appx ; Get-AppPackageLog} || $true
