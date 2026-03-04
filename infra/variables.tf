# -----------------------------------------------------------------------------
# General Variables
# -----------------------------------------------------------------------------
variable "aws_region" {
  description = "AWS 리전"
  type        = string
  default     = "ap-northeast-2"
}

variable "environment" {
  description = "환경 (dev, staging, prod)"
  type        = string
  default     = "dev"
}

variable "project_name" {
  description = "프로젝트 이름"
  type        = string
  default     = "motif-discovery"
}

# -----------------------------------------------------------------------------
# EC2 Variables
# -----------------------------------------------------------------------------
variable "instance_type" {
  description = "EC2 인스턴스 타입"
  type        = string
  default     = "r5.8xlarge"
}

variable "key_name" {
  description = "EC2 SSH 키 페어 이름"
  type        = string
  default     = null
}

variable "enable_ssh" {
  description = "SSH 접속 허용 여부"
  type        = bool
  default     = false
}

variable "allowed_ssh_cidr" {
  description = "SSH 접속 허용 CIDR (보안을 위해 본인 IP로 제한 권장)"
  type        = string
  default     = ""

  validation {
    condition     = var.allowed_ssh_cidr != "0.0.0.0/0"
    error_message = "Using 0.0.0.0/0 for allowed_ssh_cidr is not permitted. Please specify a restricted CIDR range (e.g., your IP address)."
  }
}

# -----------------------------------------------------------------------------
# S3 Variables
# -----------------------------------------------------------------------------
variable "create_s3_bucket" {
  description = "결과 저장용 S3 버킷 생성 여부"
  type        = bool
  default     = true
}

variable "s3_bucket_name" {
  description = "S3 버킷 이름 (고유해야 함)"
  type        = string
  default     = ""
}

# -----------------------------------------------------------------------------
# Spot Instance Variables
# -----------------------------------------------------------------------------
variable "use_spot_instance" {
  description = "비용 절감을 위해 스팟 인스턴스 사용 여부"
  type        = bool
  default     = false
}

variable "spot_max_price" {
  description = "스팟 인스턴스 최대 가격 (시간당 USD)"
  type        = string
  default     = "0.10"
}
