#!/bin/bash

# ==========================================
#  Git 协作助手脚本 V2.1 (兼容性与首次推送优化)
# ==========================================

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m'

get_current_branch() {
    git symbolic-ref --short HEAD 2>/dev/null || echo "DETACHED_HEAD"
}

# 智能推送函数
smart_push() {
    local branch=$1
    local force_flag=$2
    
    echo -e "${YELLOW}正在尝试推送至 origin $branch...${NC}"
    
    # 尝试推送
    if git push $force_flag origin "$branch"; then
        echo -e "${GREEN}成功推送！${NC}"
    else
        echo -e "${RED}推送失败。检测到远程可能存在本地没有的更新。${NC}"
        echo -e "${CYAN}尝试自动拉取并合并远程更改 (Allow Unrelated Histories)...${NC}"
        
        # 尝试拉取并处理可能的不相关历史
        if git pull origin "$branch" --rebase --autostash --allow-unrelated-histories; then
            echo -e "${GREEN}同步成功，再次尝试推送...${NC}"
            git push origin "$branch"
        else
            echo -e "${RED}同步失败！请手动解决冲突后再推送。${NC}"
        fi
    fi
}

# --- 主循环 ---
while true; do
    BRANCH=$(get_current_branch)
    echo -e "\n${BLUE}==============================================${NC}"
    echo -e "  Git 助手 (当前分支: ${GREEN}${BRANCH}${BLUE})  "
    echo -e "${BLUE}==============================================${NC}"
    echo "1. 查看状态"
    echo "2. 切换分支"
    echo "3. 提交并同步 (Add + Commit + Pull + Push)"
    echo "4. 仅推送"
    echo "0. 退出"
    
    # 修复某些 Shell 下 read 的兼容性
    printf "${BLUE}请输入选项: ${NC}"
    read choice

    case $choice in
        1) git status ;;
        2) 
            printf "输入目标分支名: "
            read target
            git checkout "$target" 
            ;;
        3)
            if [ -z "$(git status --porcelain)" ]; then
                echo -e "${YELLOW}没有文件需要提交。${NC}"
            else
                printf "请输入 Commit 信息: "
                read msg
                msg=${msg:-"Update"}
                git add .
                git commit -m "$msg"
            fi
            smart_push "$BRANCH"
            ;;
        4)
            smart_push "$BRANCH"
            ;;
        0) exit 0 ;;
        *) echo "无效选项" ;;
    esac

    echo -e "\n按回车键继续..."
    read dummy 
done