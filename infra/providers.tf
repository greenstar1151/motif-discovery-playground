terraform {
  required_version = ">= 1.0.0"

  required_providers {
    aws = {
      source  = "hashicorp/aws"
      version = "~> 6.30"
    }
  }

  # 선택사항: S3 백엔드 설정 (상태 파일 원격 저장)
  # backend "s3" {
  #   bucket         = "your-terraform-state-bucket"
  #   key            = "motif-discovery/terraform.tfstate"
  #   region         = "ap-northeast-2"
  #   encrypt        = true
  #   dynamodb_table = "terraform-locks"
  # }
}

provider "aws" {
  region = var.aws_region

  default_tags {
    tags = {
      Project     = "motif-discovery"
      Environment = var.environment
      ManagedBy   = "terraform"
    }
  }
}
