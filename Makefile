.PHONY: build
build:
	dotnet build src --configuration Release

.PHONY: publish
publish:
	dotnet publish src -r win10-x64 -c Release -p:PublishSingleFile=true -p:PublishTrimmed=true
	cp -f ./src/bin/Release/netcoreapp3.1/win10-x64/publish/GRAMM_Mesher.exe ../ValidationCases/Bin

.PHONY: clean
clean:
	dotnet clean src

.PHONY: test_askervein
test_askervein:
	dotnet run --project ./src "./test/Askervein_50m"
	# find "./test/Askervein_coarse" -type f -name "*.txt" | xargs rm -f
