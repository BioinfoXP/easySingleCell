#!/bin/bash

# ==========================================
#  Git 协作助手 V3.0 (速度优化版 + R包专用)
# ==========================================

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'

get_current_branch() {
    git symbolic-ref --short HEAD 2>/dev/null || echo "DETACHED_HEAD"
}

# 自动修复 R 包常见的 NAMESPACE 冲突
fix_r_namespace_conflict() {
    if [ -f "NAMESPACE" ]; then
        # 检测文件中是否有 Git 冲突标记 <<<<<<<
        if grep -q "<<<<<<<" NAMESPACE; then
            echo -e "${YELLOW}检测到 NAMESPACE 文件存在冲突，正在自动重置...${NC}"
            # 策略：直接删除 NAMESPACE，让 R (devtools::document) 重新生成
            rm NAMESPACE
            echo -e "${GREEN}NAMESPACE 已删除，请在 R 中运行 document() 重新生成。${NC}"
        fi
    fi
}

# 优化的同步函数：先拉后推
smart_sync() {
    local branch=$1
    echo -e "${CYAN}--- 步骤 1: 从远程拉取更新 (Rebase) ---${NC}"
    
    # 使用 --rebase 保持提交历史整洁，--autostash 暂存未提交的修改
    if git pull origin "$branch" --rebase --autostash; then
        echo -e "${GREEN}拉取成功 (或已是最新)${NC}"
    else
        echo -e "${RED}拉取遇到冲突！${NC}"
        fix_r_namespace_conflict
        echo -e "${YELLOW}尝试解决冲突后，请手动运行 git rebase --continue${NC}"
        return 1
    fi

    echo -e "${CYAN}--- 步骤 2: 推送到远程 ---${NC}"
    if git push origin "$branch"; then
        echo -e "${GREEN}>>> 推送成功！ <<<${NC}"
    else
        echo -e "${RED}推送失败，请检查网络或权限。${NC}"
    fi
}

# --- 主循环 ---
while true; do
    BRANCH=$(get_current_branch)
    echo -e "\n${BLUE}==============================================${NC}"
    echo -e "  Git 极速助手 (当前分支: ${GREEN}${BRANCH}${BLUE})  "
    echo -e "${BLUE}==============================================${NC}"
    echo "1. 查看状态 (git status)"
    echo "2. 切换分支"
    echo "3. 提交并同步 (推荐: Add+Commit -> Pull -> Push)"
    echo "4. 强制推送 (慎用)"
    echo "0. 退出"
    
    printf "${BLUE}请输入选项: ${NC}"
    read -r choice  # 使用 -r 防止反斜杠转义

    case $choice in
        1) git status ;;
        2) 
            printf "输入目标分支名: "
            read -r target
            git checkout "$target" 
            ;;
        3)
            # 1. 检查是否有文件需要 Add
            if [ -n "$(git status --porcelain)" ]; then
                printf "请输入 Commit 信息 (回车默认 'Update'): "
                read -r msg
                msg=${msg:-"Update"}
                
                # 特殊处理：如果是 R 包开发，建议先检查 NAMESPACE
                fix_r_namespace_conflict
                
                git add .
                git commit -m "$msg"
            else
                echo -e "${YELLOW}工作区干净，无新文件提交。直接进入同步步骤...${NC}"
            fi
            
            # 2. 执行同步 (先拉后推)
            smart_sync "$BRANCH"
            ;;
        4)
            echo -e "${RED}警告: 强制推送可能会覆盖远程他人的代码！${NC}"
            printf "确认强制推送? (y/n): "
            read -r confirm
            if [ "$confirm" = "y" ]; then
                git push origin "$BRANCH" --force
            fi
            ;;
        0) exit 0 ;;
        *) echo "无效选项" ;;
    esac

    echo -e "\n按回车键继续..."
    read -r dummy 
done