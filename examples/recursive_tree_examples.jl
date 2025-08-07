# 递归构造二叉树的完整示例
# 这个文件展示了在Julia中使用递归构造二叉树的各种方法

using TreeContractor

println("=== 递归构造二叉树示例 ===\n")

# 示例1: 从排序数组构造平衡二叉搜索树
println("1. 从排序数组构造平衡二叉搜索树")
println("输入: [1, 2, 3, 4, 5, 6, 7]")
sorted_arr = [1, 2, 3, 4, 5, 6, 7]
tree1 = construct_bst_from_sorted_array(sorted_arr)

println("构造的树结构:")
println("     4")
println("    / \\")
println("   2   6")
println("  / \\ / \\")
println(" 1  3 5  7")

println("中序遍历结果: ", inorder_traversal(tree1))
println("先序遍历结果: ", preorder_traversal(tree1))
println("后序遍历结果: ", postorder_traversal(tree1))
println("树的高度: ", tree_height_recursive(tree1))
println("树的节点数: ", tree_size_recursive(tree1))
println()

# 示例2: 从先序遍历和中序遍历构造二叉树
println("2. 从先序遍历和中序遍历构造二叉树")
println("先序遍历: [3, 9, 20, 15, 7]")
println("中序遍历: [9, 3, 15, 20, 7]")

preorder = [3, 9, 20, 15, 7]
inorder = [9, 3, 15, 20, 7]
tree2 = construct_from_preorder_inorder(preorder, inorder)

println("构造的树结构:")
println("     3")
println("    / \\")
println("   9  20")
println("      / \\")
println("     15  7")

println("验证先序遍历: ", preorder_traversal(tree2))
println("验证中序遍历: ", inorder_traversal(tree2))
println()

# 示例3: 从后序遍历和中序遍历构造二叉树
println("3. 从后序遍历和中序遍历构造二叉树")
println("后序遍历: [9, 15, 7, 20, 3]")
println("中序遍历: [9, 3, 15, 20, 7]")

postorder = [9, 15, 7, 20, 3]
inorder2 = [9, 3, 15, 20, 7]
tree3 = construct_from_postorder_inorder(postorder, inorder2)

println("构造的树结构:")
println("     3")
println("    / \\")
println("   9  20")
println("      / \\")
println("     15  7")

println("验证后序遍历: ", postorder_traversal(tree3))
println("验证中序遍历: ", inorder_traversal(tree3))
println()

# 示例4: 从层序遍历构造二叉树
println("4. 从层序遍历构造二叉树")
println("层序遍历: [3, 9, 20, nothing, nothing, 15, 7]")

level_order = [3, 9, 20, nothing, nothing, 15, 7]
tree4 = construct_from_level_order(level_order)

println("构造的树结构:")
println("     3")
println("    / \\")
println("   9  20")
println("      / \\")
println("     15  7")

println("层序遍历结果: ", level_order_traversal(tree4))
println()

# 示例5: 递归插入构造二叉搜索树
println("5. 递归插入构造二叉搜索树")
println("插入顺序: [5, 3, 7, 1, 4, 6, 8]")

tree5 = nothing
values = [5, 3, 7, 1, 4, 6, 8]

for val in values
    tree5 = insert_bst_recursive(tree5, val)
end

println("构造的树结构:")
println("     5")
println("    / \\")
println("   3   7")
println("  / \\ / \\")
println(" 1  4 6  8")

println("中序遍历结果: ", inorder_traversal(tree5))
println("验证是否为有效BST: ", is_valid_bst(tree5))
println()

# 示例6: 递归删除节点
println("6. 递归删除节点")
println("原始树:")
println("     5")
println("    / \\")
println("   3   7")
println("  / \\ / \\")
println(" 1  4 6  8")

# 删除叶子节点
tree6 = deepcopy(tree5)
tree6 = delete_bst_recursive(tree6, 1)
println("删除节点1后:")
println("     5")
println("    / \\")
println("   3   7")
println("    \\ / \\")
println("    4 6  8")

# 删除只有一个子节点的节点
tree6 = delete_bst_recursive(tree6, 3)
println("删除节点3后:")
println("     5")
println("    / \\")
println("   4   7")
println("      / \\")
println("     6  8")

# 删除有两个子节点的节点
tree6 = delete_bst_recursive(tree6, 5)
println("删除节点5后:")
println("     6")
println("    / \\")
println("   4   7")
println("        \\")
println("         8")

println("最终中序遍历: ", inorder_traversal(tree6))
println()

# 示例7: 递归查找节点
println("7. 递归查找节点")
tree7 = deepcopy(tree5)
println("在树中查找节点4: ", find_node_recursive(tree7, 4) !== nothing ? "找到" : "未找到")
println("在树中查找节点10: ", find_node_recursive(tree7, 10) !== nothing ? "找到" : "未找到")
println()

# 示例8: 性能比较
println("8. 不同构造方法的性能比较")
using BenchmarkTools

# 测试从排序数组构造的性能
sorted_test = collect(1:1000)
println("从排序数组构造1000个节点的树:")
@btime construct_bst_from_sorted_array($sorted_test)

# 测试递归插入的性能
println("递归插入1000个随机节点:")
random_test = rand(1:1000, 1000)
@btime begin
    tree = nothing
    for val in random_test
        tree = insert_bst_recursive(tree, val)
    end
end

println("\n=== 递归构造二叉树的总结 ===")
println("1. 从排序数组构造: 时间复杂度O(n), 空间复杂度O(log n)")
println("2. 从遍历序列构造: 时间复杂度O(n), 空间复杂度O(n)")
println("3. 递归插入构造: 时间复杂度O(n log n), 空间复杂度O(n)")
println("4. 递归删除: 时间复杂度O(h), 其中h是树的高度")
println("5. 递归查找: 时间复杂度O(h), 其中h是树的高度")
println("6. 递归验证: 时间复杂度O(n), 空间复杂度O(h)") 