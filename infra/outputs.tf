# -----------------------------------------------------------------------------
# Outputs
# -----------------------------------------------------------------------------
output "vpc_id" {
  description = "VPC ID"
  value       = aws_vpc.main.id
}

output "instance_id" {
  description = "EC2 인스턴스 ID"
  value       = aws_instance.motif_discovery.id
}

output "instance_public_ip" {
  description = "EC2 인스턴스 퍼블릭 IP"
  value       = aws_instance.motif_discovery.public_ip
}

output "instance_public_dns" {
  description = "EC2 인스턴스 퍼블릭 DNS"
  value       = aws_instance.motif_discovery.public_dns
}

output "ssh_command" {
  description = "SSH 접속 명령어"
  value       = var.key_name != null ? "ssh -i ~/.ssh/${var.key_name}.pem ubuntu@${aws_instance.motif_discovery.public_ip}" : "SSH key not configured. Use SSM Session Manager instead."
}

output "ssm_command" {
  description = "SSM 세션 매니저 접속 명령어"
  value       = "aws ssm start-session --target ${aws_instance.motif_discovery.id}"
}

output "s3_bucket_name" {
  description = "결과 저장용 S3 버킷 이름"
  value       = var.create_s3_bucket ? aws_s3_bucket.results[0].id : "S3 bucket not created"
}

output "s3_bucket_arn" {
  description = "S3 버킷 ARN"
  value       = var.create_s3_bucket ? aws_s3_bucket.results[0].arn : "S3 bucket not created"
}

# 스팟 인스턴스 정보 (사용 시)
output "spot_instance_id" {
  description = "스팟 인스턴스 ID"
  value       = var.use_spot_instance ? aws_spot_instance_request.motif_discovery_spot[0].spot_instance_id : "Spot instance not requested"
}

output "spot_instance_public_ip" {
  description = "스팟 인스턴스 퍼블릭 IP"
  value       = var.use_spot_instance ? aws_spot_instance_request.motif_discovery_spot[0].public_ip : "Spot instance not requested"
}
