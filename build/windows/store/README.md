# GIMP MSIX HowTo

## Stable and Development releases

The MS Store/Partner Center hosts our stable and development point releases:
https://partner.microsoft.com/dashboard/apps-and-games/overview

Base rule to update the "GIMP" or "GIMP (Preview)" entry:

* Regularly, a .msixupload will be generated at each tagged commit. It also can be
  with a custom pipeline valuating `GIMP_CI_MS_STORE` variable to `${REVISION}` as
  briefly explained in ["Versioning the MSIX"](#versioning-the-msix) section.
  In the process, it will be auto submitted (without changelog) to Partner Center.

## Maintaining the MSIX

If the submission from CI does not work, check on MS Entra Admin Center for the
client secret of 'GITLAB_CI' app. If expired, create a new one and add it on
GitLab as a masked, hidden and protected variable. Then, retrigger the Store job.
That is usually faster than doing a manual submission as per instructions below:

1. The link above will open 'Partner Center', the frontend to submit .msixupload

2. With your Microsoft Entra ID account (@*onmicrosoft.com), log-in as customary

3. Clicking on the "Start update" button will make possible to change some things.
   Only 'Packages' and 'Store listings' sections are needed. On 'Packages' you will
   add the generated .msixupload and on 'Store listings' the brief changelog.

If the .msix* starts to be refused to certification or to self-signing,
run `build\windows\store\3_dist-gimp-winsdk.ps1 WACK` locally to see if it
still complies with the latest Windows policies. Make sure to update WinSDK.

If the .msix* starts to be refused to self-signing due to the .pfx file, then
generate a new one with the commands below and commit it to this dir.

```pwsh
$pseudo_gimp = "pseudo-gimp_$(Get-Date -Format yyyy)"
```

```pwsh
New-SelfSignedCertificate -Type Custom -Subject "$(([xml](Get-Content build\windows\store\AppxManifest.xml)).Package.Identity.Publisher)" -KeyUsage DigitalSignature -FriendlyName "$pseudo_gimp" -CertStoreLocation "Cert:\CurrentUser\My" -TextExtension @("2.5.29.37={text}1.3.6.1.5.5.7.3.3", "2.5.29.19={text}")
```

```pwsh
Export-PfxCertificate -Cert "Cert:\CurrentUser\My\$(Get-ChildItem Cert:\CurrentUser\My | Where-Object FriendlyName -EQ "$pseudo_gimp" | Select-Object -ExpandProperty Thumbprint)" -FilePath "${pseudo_gimp}.pfx" -Password (ConvertTo-SecureString -String eek -Force -AsPlainText)
```

## Versioning the MSIX

* Every new .msixupload submission (with different content) needs a bumped version.
  But MS Store reserves revisions (the fourth digit after last dot) for its own use.
  https://learn.microsoft.com/en-us/windows/apps/publish/publish-your-app/msix/app-package-requirements
  So, our MSIX have a special micro versioning to we have revisions. Check the .ps1.
