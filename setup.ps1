if (Test-Path build) {
    Remove-Item -Recurse -Force build
}

New-Item -ItemType Directory -Force -Path build
Set-Location build

cmake -G "Ninja" -DCMAKE_BUILD_TYPE=Release ..

ninja

Copy-Item bin/Box_test.exe ../out/

cd ../out/




