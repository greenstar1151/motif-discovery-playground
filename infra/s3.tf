# -----------------------------------------------------------------------------
# S3 Bucket for Results
# -----------------------------------------------------------------------------
resource "random_id" "bucket_suffix" {
  count       = var.create_s3_bucket && var.s3_bucket_name == "" ? 1 : 0
  byte_length = 4
}

resource "aws_s3_bucket" "results" {
  count = var.create_s3_bucket ? 1 : 0

  bucket = var.s3_bucket_name != "" ? var.s3_bucket_name : "${var.project_name}-results-${random_id.bucket_suffix[0].hex}"

  tags = {
    Name = "${var.project_name}-results"
  }
}

resource "aws_s3_bucket_versioning" "results" {
  count  = var.create_s3_bucket ? 1 : 0
  bucket = aws_s3_bucket.results[0].id

  versioning_configuration {
    status = "Enabled"
  }
}

resource "aws_s3_bucket_server_side_encryption_configuration" "results" {
  count  = var.create_s3_bucket ? 1 : 0
  bucket = aws_s3_bucket.results[0].id

  rule {
    apply_server_side_encryption_by_default {
      sse_algorithm = "AES256"
    }
  }
}

resource "aws_s3_bucket_public_access_block" "results" {
  count  = var.create_s3_bucket ? 1 : 0
  bucket = aws_s3_bucket.results[0].id

  block_public_acls       = true
  block_public_policy     = true
  ignore_public_acls      = true
  restrict_public_buckets = true
}

# 오래된 결과 자동 삭제를 위한 라이프사이클 규칙
resource "aws_s3_bucket_lifecycle_configuration" "results" {
  count  = var.create_s3_bucket ? 1 : 0
  bucket = aws_s3_bucket.results[0].id

  rule {
    id     = "cleanup-old-results"
    status = "Enabled"

    # 90일 후 오래된 버전 삭제
    noncurrent_version_expiration {
      noncurrent_days = 90
    }

    # 불완전한 멀티파트 업로드 정리
    abort_incomplete_multipart_upload {
      days_after_initiation = 7
    }
  }
}
