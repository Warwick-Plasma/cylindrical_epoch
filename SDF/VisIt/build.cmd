xml2info.exe  -clobber SDF2.xml
xml2cmake.exe -clobber SDF2.xml

"C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe" -G "Visual Studio 15 2017 Win64" .

"C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\Common7\IDE\devenv.exe" /Build Release SDF_database.sln
