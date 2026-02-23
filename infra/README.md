# Motif Discovery - AWS Infrastructure

`infra/`는 Terraform으로 Motif Discovery 실행용 AWS 인프라를 생성합니다.

## 현재 구성

- VPC + Public Subnet + Internet Gateway + Route Table
- Security Group (옵션 SSH + 전체 egress)
- EC2 (Ubuntu 22.04 AMI, 기본 30GB gp3)
- IAM Role/Instance Profile (SSM 기본 포함)
- S3 버킷 (선택: 버전닝, SSE, 퍼블릭 차단, 라이프사이클)
- Spot Instance Request (선택)

> 참고: 현재 코드 기준으로는 일반 EC2 인스턴스가 항상 생성되고,
> `use_spot_instance = true`일 때 Spot 인스턴스가 추가로 생성됩니다.

## 사전 요구사항

- Terraform >= 1.0
- AWS CLI 및 자격 증명 설정

```bash
aws configure
```

## 빠른 시작

```bash
cd infra
cp terraform.tfvars.example terraform.tfvars
# 필요 값 수정 (예: key_name, allowed_ssh_cidr, instance_type)

terraform init
terraform plan
terraform apply
```

`prod.tfvars`를 사용하려면:

```bash
terraform plan -var-file=prod.tfvars
terraform apply -var-file=prod.tfvars
```

## 접속

```bash
# 출력값 확인
terraform output

# SSH (key_name 사용 시)
ssh -i ~/.ssh/<key_name>.pem ubuntu@<instance_public_ip>

# SSM
aws ssm start-session --target <instance_id>
```

## 인스턴스 준비 상태

`scripts/user_data.sh`에서 아래 작업을 자동 수행합니다.

- 필수 패키지 설치 (cmake, git, awscli 등)
- 저장소 클론: `~/motif-discovery-playground`

이후 빌드/실행은 인스턴스에서 수동으로 진행합니다.

```bash
cd ~/motif-discovery-playground
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
./build/motif_benchmark --large
```

## 주요 변수

- `instance_type` (기본: `r5.8xlarge`)
- `key_name`
- `enable_ssh`, `allowed_ssh_cidr`
- `create_s3_bucket`, `s3_bucket_name`
- `use_spot_instance`, `spot_max_price`

## 정리

```bash
terraform destroy
```
