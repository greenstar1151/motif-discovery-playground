#!/bin/bash
set -e

# 로그 설정
exec > >(tee /var/log/user-data.log) 2>&1
echo "=== Starting user data script at $(date) ==="

# -----------------------------------------------------------------------------
# 시스템 패키지 업데이트 및 필수 도구 설치
# -----------------------------------------------------------------------------
echo "=== Installing system packages ==="
apt-get update
apt-get install -y \
    build-essential \
    cmake \
    git \
    htop \
    awscli \
    unzip \
    jq \
    python3 \
    python3-pip

# -----------------------------------------------------------------------------
# 프로젝트 클론
# -----------------------------------------------------------------------------
echo "=== Cloning project ==="
WORK_DIR="/home/ubuntu/motif-discovery-playground"

# GitHub에서 프로젝트 클론
git clone https://github.com/greenstar1151/motif-discovery-playground.git "$WORK_DIR"

# 소유권 설정
chown -R ubuntu:ubuntu "$WORK_DIR"

# -----------------------------------------------------------------------------
# 환경 변수 설정
# -----------------------------------------------------------------------------
echo "=== Setting up environment ==="
cat >> /home/ubuntu/.bashrc << 'EOF'

# Motif Discovery 환경 변수
export MOTIF_HOME="/home/ubuntu/motif-discovery-playground"
export PATH="$MOTIF_HOME/build:$PATH"

# 프로젝트 디렉토리로 이동
cd "$MOTIF_HOME"
EOF

# -----------------------------------------------------------------------------
# 완료 메시지
# -----------------------------------------------------------------------------
echo "=== User data script completed at $(date) ==="
echo "Project cloned at: $WORK_DIR"
echo ""
echo "Next steps after SSH login:"
echo "  1. cd ~/motif-discovery-playground"
echo "  2. cmake -B build -DCMAKE_BUILD_TYPE=Release"
echo "  3. cmake --build build -j\$(nproc)"
echo "  4. Prepare your data"
echo "  5. Run your experiments"
