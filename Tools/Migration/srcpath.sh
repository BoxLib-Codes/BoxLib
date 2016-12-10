#!/usr/bin/env bash

OLD_SRC_PATH="Src\/C_AMRLib"
NEW_SRC_PATH="Src\/Amr"
find . -path .git -prune -o -type f -exec grep -Iq . {} \; -exec sed -i 's/'"${OLD_SRC_PATH}"'/'"${NEW_SRC_PATH}"'/g' {} +

OLD_SRC_PATH="Src\/C_AmrCoreLib"
NEW_SRC_PATH="Src\/AmrCore"
find . -path .git -prune -o -type f -exec grep -Iq . {} \; -exec sed -i 's/'"${OLD_SRC_PATH}"'/'"${NEW_SRC_PATH}"'/g' {} +

OLD_SRC_PATH="Src\/C_BaseLib"
NEW_SRC_PATH="Src\/Base"
find . -path .git -prune -o -type f -exec grep -Iq . {} \; -exec sed -i 's/'"${OLD_SRC_PATH}"'/'"${NEW_SRC_PATH}"'/g' {} +

OLD_SRC_PATH="Src\/C_BoundaryLib"
NEW_SRC_PATH="Src\/Boundary"
find . -path .git -prune -o -type f -exec grep -Iq . {} \; -exec sed -i 's/'"${OLD_SRC_PATH}"'/'"${NEW_SRC_PATH}"'/g' {} +

OLD_SRC_PATH="Src\/C_ParticleLib"
NEW_SRC_PATH="Src\/Particle"
find . -path .git -prune -o -type f -exec grep -Iq . {} \; -exec sed -i 's/'"${OLD_SRC_PATH}"'/'"${NEW_SRC_PATH}"'/g' {} +

