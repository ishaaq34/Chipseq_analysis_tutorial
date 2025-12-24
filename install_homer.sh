#!/bin/bash

# HOMER Installation Script for Mac
# Navigate to your project directory first

# Create tools directory
mkdir -p tools
cd tools

# Download HOMER configuration script
curl -O http://homer.ucsd.edu/homer/configureHomer.pl

# Download HOMER v4.11 package
curl -L http://homer.ucsd.edu/homer/data/software/homer.v4.11.zip -o homer.zip

# Extract HOMER
unzip homer.zip

# Remove zip file
rm homer.zip

# Install/configure HOMER (compiles C++ components)
perl configureHomer.pl -install

# Go back to project root
cd ..

# Add HOMER to PATH permanently
HOMER_PATH="export PATH=\$PATH:$(pwd)/tools/bin/"
echo ""
echo "Adding HOMER to your PATH in ~/.zshrc..."

# Check if already in zshrc
if grep -q "tools/bin/" ~/.zshrc 2>/dev/null; then
    echo "HOMER PATH already exists in ~/.zshrc"
else
    echo "$HOMER_PATH" >> ~/.zshrc
    echo "✅ Added HOMER to ~/.zshrc"
fi

# Load HOMER for current session
export PATH=$PATH:$(pwd)/tools/bin/

# Test installation
echo ""
echo "Testing HOMER installation..."
findMotifsGenome.pl 2>&1 | head -5

echo ""
echo "✅ HOMER installed successfully!"
echo ""
echo "HOMER is now available globally. You can:"
echo "  1. Reload your shell: source ~/.zshrc"
echo "  2. Or open a new terminal window"
echo ""
echo "Test it by running: findMotifsGenome.pl"
