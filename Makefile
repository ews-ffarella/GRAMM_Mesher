.PHONY: build
build:
	dotnet build src --configuration Release

.PHONY: publish
publish:
	dotnet publish src -r win10-x64 -c Release -p:PublishSingleFile=true -p:PublishTrimmed=true

.PHONY: clean
clean:
	dotnet clean src

.PHONY: test_a_coarse
test_a_coarse:
	dotnet run --project ./src "./test/Askervein_coarse"
	# find "./test/Askervein_coarse" -type f -name "*.txt" | xargs rm -f


.PHONY: test_a_fine
test_a_fine:
	dotnet run --project ./src "./test/Askervein_fine"
	# find "./test/Askervein_coarse" -type f -name "*.txt" | xargs rm -f