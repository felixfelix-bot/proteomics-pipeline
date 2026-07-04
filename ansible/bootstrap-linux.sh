#!/usr/bin/env bash
###############################################################################
# bootstrap-linux.sh
# Installs Ansible on the Linux control machine so it can provision targets.
# Run this ONCE on the Linux machine that will control provisioning.
###############################################################################
set -euo pipefail

echo "=========================================="
echo " Proteomics Pipeline — Ansible Bootstrap"
echo " Linux Control Machine Setup"
echo "=========================================="

# Detect distro
if [ -f /etc/debian_version ]; then
    DISTRO="debian"
    echo "Detected Debian/Ubuntu..."
elif [ -f /etc/fedora-release ]; then
    DISTRO="fedora"
    echo "Detected Fedora/RHEL..."
else
    echo "Unsupported distro. Install Ansible manually:"
    echo "  pip3 install ansible"
    exit 1
fi

# Install Python + pip if missing
if ! command -v python3 &>/dev/null; then
    echo "Installing Python3..."
    if [ "$DISTRO" = "debian" ]; then
        sudo apt-get update && sudo apt-get install -y python3 python3-pip python3-venv
    else
        sudo dnf install -y python3 python3-pip
    fi
fi

# Install Ansible
if ! command -v ansible &>/dev/null; then
    echo "Installing Ansible..."
    if [ "$DISTRO" = "debian" ]; then
        sudo apt-get update
        sudo apt-get install -y ansible
    else
        sudo dnf install -y ansible
    fi
fi

# Install pywinrm for Windows targets
echo "Installing pywinrm (for Windows target management)..."
pip3 install --user pywinrm 2>/dev/null || sudo pip3 install pywinrm

# Verify
echo ""
echo "=========================================="
echo " Verification:"
echo "=========================================="
ansible --version | head -3
echo ""

echo "Ansible is ready. Next steps:"
echo ""
echo "1. Edit inventory.yml to set your Windows machine IP/user/password:"
echo "   vim ansible/inventory.yml"
echo ""
echo "2. On the Windows machine, run bootstrap-windows.ps1 as Administrator"
echo "   to enable WinRM."
echo ""
echo "3. Run the playbook:"
echo "   cd ansible && ansible-playbook -i inventory.yml setup.yml"
echo ""
echo "To provision THIS Linux machine locally:"
echo "   cd ansible && ansible-playbook -i inventory.yml setup.yml -l linux_targets"
