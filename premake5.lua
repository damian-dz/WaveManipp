workspace "WaveManipp"
    configurations { "DynDebug", "DynRelease", "StaDebug", "StaRelease" }
    platforms { "x86", "x64" }
    
    filter "platforms:x86"
        architecture "x86"

    filter "platforms:x64"
        architecture "x86_64"
        
    filter "configurations:Dyn*"
        kind "SharedLib"
        defines "WAVEMANIPP_DLL"
        
    filter "configurations:Sta*"
        kind "StaticLib"

project "WaveManipp"
    language "C++"
    location "WaveManipp"
    targetdir "bin/%{cfg.platform}/%{cfg.buildcfg}"
    files { "api/**.hpp", "src/**.cpp" }
    includedirs "api/"

    filter "system:windows"       
        cppdialect "C++17"
        defines "_CRT_SECURE_NO_WARNINGS"
        systemversion "latest"
        filter "configurations:Dyn*"
            staticruntime "Off"
        filter "configurations:Sta*"
            staticruntime "On"   

    filter "configurations:*Debug"
        defines "DEBUG"
        symbols "On"

    filter "configurations:*Release"
        defines "NDEBUG"
        optimize "On"
