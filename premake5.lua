workspace "WaveManipp"
    architecture "x64"
    configurations { "Debug", "Release" }

project "WaveManipp"
    language "C++"
	location "WaveManipp"
    kind "StaticLib"
	targetdir "bin/%{cfg.buildcfg}"
	files {
        "src/**.h",
        "src/**.hpp",
        "src/**.cpp"
    }
	
	defines { }

	filter "system:windows"       
	    cppdialect "C++17"
		defines {
		    "_CRT_SECURE_NO_WARNINGS"
		} 
		staticruntime "On"
		systemversion "latest"

    filter "configurations:Debug"
	    defines { "DEBUG" }
        symbols "On"

    filter "configurations:Release"
        defines { "NDEBUG" }
        optimize "On"