# -----------------------------------------------------------------------------
# Data Sources
# -----------------------------------------------------------------------------

# Amazon Linux 2023 AMI (x86_64)
data "aws_ami" "amazon_linux_2023" {
  most_recent = true
  owners      = ["amazon"]

  filter {
    name   = "name"
    values = ["al2023-ami-*-x86_64"]
  }

  filter {
    name   = "virtualization-type"
    values = ["hvm"]
  }
}

# Ubuntu 22.04 LTS AMI (대안)
data "aws_ami" "ubuntu" {
  most_recent = true
  owners      = ["099720109477"] # Canonical

  filter {
    name   = "name"
    values = ["ubuntu/images/hvm-ssd/ubuntu-jammy-22.04-amd64-server-*"]
  }

  filter {
    name   = "virtualization-type"
    values = ["hvm"]
  }
}

# -----------------------------------------------------------------------------
# IAM Role for EC2
# -----------------------------------------------------------------------------
resource "aws_iam_role" "ec2_role" {
  name = "${var.project_name}-ec2-role"

  assume_role_policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Action = "sts:AssumeRole"
        Effect = "Allow"
        Principal = {
          Service = "ec2.amazonaws.com"
        }
      }
    ]
  })
}

# S3 접근 정책 (결과 업로드용)
resource "aws_iam_role_policy" "s3_access" {
  count = var.create_s3_bucket ? 1 : 0
  name  = "${var.project_name}-s3-access"
  role  = aws_iam_role.ec2_role.id

  policy = jsonencode({
    Version = "2012-10-17"
    Statement = [
      {
        Effect = "Allow"
        Action = [
          "s3:PutObject",
          "s3:GetObject",
          "s3:ListBucket"
        ]
        Resource = [
          aws_s3_bucket.results[0].arn,
          "${aws_s3_bucket.results[0].arn}/*"
        ]
      }
    ]
  })
}

# SSM 세션 매니저 접근 (SSH 없이 접속 가능)
resource "aws_iam_role_policy_attachment" "ssm" {
  role       = aws_iam_role.ec2_role.name
  policy_arn = "arn:aws:iam::aws:policy/AmazonSSMManagedInstanceCore"
}

resource "aws_iam_instance_profile" "ec2_profile" {
  name = "${var.project_name}-ec2-profile"
  role = aws_iam_role.ec2_role.name
}

# -----------------------------------------------------------------------------
# EC2 Instance
# -----------------------------------------------------------------------------
resource "aws_instance" "motif_discovery" {
  ami                    = data.aws_ami.ubuntu.id
  instance_type          = var.instance_type
  key_name               = var.key_name
  subnet_id              = aws_subnet.public.id
  vpc_security_group_ids = [aws_security_group.motif_discovery.id]
  iam_instance_profile   = aws_iam_instance_profile.ec2_profile.name

  root_block_device {
    volume_size           = 30
    volume_type           = "gp3"
    delete_on_termination = true
    encrypted             = true
  }

  user_data = file("${path.module}/scripts/user_data.sh")

  tags = {
    Name = "${var.project_name}-instance"
  }
}

# -----------------------------------------------------------------------------
# Spot Instance (비용 절감 옵션)
# -----------------------------------------------------------------------------
resource "aws_spot_instance_request" "motif_discovery_spot" {
  count = var.use_spot_instance ? 1 : 0

  ami                    = data.aws_ami.ubuntu.id
  instance_type          = var.instance_type
  key_name               = var.key_name
  subnet_id              = aws_subnet.public.id
  vpc_security_group_ids = [aws_security_group.motif_discovery.id]
  iam_instance_profile   = aws_iam_instance_profile.ec2_profile.name

  spot_price           = var.spot_max_price
  wait_for_fulfillment = true
  spot_type            = "one-time"

  root_block_device {
    volume_size           = 50
    volume_type           = "gp3"
    delete_on_termination = true
    encrypted             = true
  }

  user_data = templatefile("${path.module}/scripts/user_data.sh", {
    s3_bucket = var.create_s3_bucket ? aws_s3_bucket.results[0].id : ""
  })

  tags = {
    Name = "${var.project_name}-spot-instance"
  }
}
