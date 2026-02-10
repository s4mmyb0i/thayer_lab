# Data Structures and Algorithms Tutorial

Over the course of the semester, I worked through the NeetCode 150 problem set as a structured way to strengthen my understanding of core data structures and algorithms, particularly for software engineering interview preparation. While the preexisting Algorithms course focuses heavily on theoretical analysis, it does not emphasize the practical use and implementation of common data structures. NeetCode provided a comprehensive and systematic framework that helped bridge that gap and gave me repeated hands-on practice applying these concepts.

Arrays and hashing was one of the areas where I saw the most improvement. Early on, I had conceptual gaps in how to effectively use hash maps and frequency counting, but working through problems like Majority Element and Longest Common Prefix made these techniques feel much more natural. I became more comfortable identifying when constant-time lookups could simplify a problem and improve time complexity. This process also helped me better understand Python-specific details, such as the fact that ```list.pop()``` is an $O(1)$ operation only when removing from the end of a list, but becomes $O(n)$ when removing from other positions—an important consideration when writing efficient code.

A major learning curve for me was two-pointer techniques. Initially, these were challenging, but through repeated exposure I learned several common patterns, including:

- Opposite-direction pointers
- Same-direction sliding pointers
- Fast and slow pointers

Understanding when and why to use each of these patterns significantly improved my ability to reason about array and string problems efficiently.

Although I had previously studied search algorithms such as binary search in my algorithms course, NeetCode helped me become fluent in implementing them correctly and efficiently. Problems like Binary Search, Search in Rotated Sorted Array, and Find Minimum in Rotated Sorted Array strengthened my intuition for maintaining invariants and handling edge cases.

I also developed stronger skills in one-dimensional dynamic programming, particularly through problems like Climbing Stairs, House Robber, and Maximum Subarray. These problems helped me learn how to define state, identify recurrence relations, and optimize space usage when possible. I have also begun working on two-dimensional dynamic programming, which has further expanded my problem-solving toolkit.

Linked lists were initially difficult, especially when it came to traversal and recognizing when they are the appropriate data structure for a given problem. However, problems such as Reverse Linked List and Merge Two Sorted Lists made pointer manipulation clearer and increased my confidence in working with linked structures.

Additionally, I became more comfortable with greedy algorithms and stack-based problems. Examples like Valid Parentheses, Min Stack, and Daily Temperatures demonstrated how stacks can efficiently track state, while greedy problems reinforced the idea that locally optimal decisions can lead to globally correct solutions.

Overall, completing NeetCode 150 gave me a strong foundational structure for approaching data structure and algorithm problems. It complemented my previous coursework by emphasizing implementation, pattern recognition, and problem-solving under constraints—skills that are essential for both software engineering interviews and practical development work.