# #' Generate Word Cloud from Text Sentences
# #'
# #' @param sentences Character vector of input sentences
# #' @param max_words Maximum number of words to display (default 100)
# #' @param remove_stopwords Logical, whether to remove stopwords (default TRUE)
# #' @param stopwords Custom stopwords vector (default NULL uses common English stopwords)
# #' @param random_order Plot words in random order (default FALSE)
# #' @param min_freq Minimum word frequency to include (default 1)
# #' @param colors Color palette for words (default rainbow palette)
# #' @param bg_color Background color (default "white")
# #' @param ... Additional arguments passed to ggplot2::theme()
# #'
# #' @return ggplot object containing the word cloud
# wordcloud_generator <- function(sentences, 
#                                max_words = 100,
#                                remove_stopwords = TRUE,
#                                stopwords = NULL,
#                                random_order = FALSE,
#                                min_freq = 1,
#                                colors = grDevices::rainbow(max_words),
#                                bg_color = "white",
#                                ...) {
  
#   # Basic text processing function
#   process_text <- function(text) {
#     text <- tolower(text)
#     text <- gsub("[[:punct:]]", "", text)
#     text <- gsub("[[:digit:]]", "", text)
#     text <- gsub("\\s+", " ", text)
#     trimws(text)
#   }
  
#   # Set default stopwords if not provided
#   if (is.null(stopwords)) {
#     stopwords <- c("a", "an", "the", "and", "or", "but", "to", "of", "for", 
#                   "in", "on", "at", "by", "with", "as", "is", "are", "was", 
#                   "were", "be", "been", "being", "have", "has", "had", "do", 
#                   "does", "did", "will", "would", "should", "can", "could", 
#                   "may", "might", "must", "i", "you", "he", "she", "it", 
#                   "we", "they", "me", "him", "her", "us", "them", "this", 
#                   "that", "these", "those", "my", "your", "his", "its", 
#                   "our", "their")
#   }
  
#   # Process all input sentences
#   processed_text <- process_text(paste(sentences, collapse = " "))
#   words <- unlist(strsplit(processed_text, " "))
  
#   # Filter words
#   if (remove_stopwords) {
#     words <- words[!words %in% stopwords]
#   }
#   words <- words[nchar(words) > 0]
  
#   # Calculate word frequencies
#   word_freq <- sort(table(words), decreasing = TRUE)
#   word_freq <- word_freq[word_freq >= min_freq]
#   if (length(word_freq) > max_words) {
#     word_freq <- word_freq[1:max_words]
#   }
  
#   # Prepare plot data
#   plot_data <- data.frame(
#     word = names(word_freq),
#     freq = as.numeric(word_freq),
#     stringsAsFactors = FALSE
#   )
  
#   # Word layout function
#   layout_words <- function(words, freqs) {
#     n <- length(words)
#     sizes <- sqrt(freqs/max(freqs)) * 10
    
#     theta <- seq(0, 2*pi, length.out = n + 1)[1:n]
#     if (random_order) theta <- sample(theta)
    
#     r <- seq(1, 5, length.out = n)
#     data.frame(
#       word = words,
#       x = r * cos(theta),
#       y = r * sin(theta),
#       size = sizes,
#       angle = runif(n, -30, 30)
#     )
#   }
  
#   # Generate layout
#   layout <- layout_words(plot_data$word, plot_data$freq)
  
#   # Create ggplot
#   p <- ggplot2::ggplot(layout, ggplot2::aes(x = x, y = y)) +
#     ggplot2::geom_text(
#       ggplot2::aes(
#         label = word,
#         size = size,
#         angle = angle,
#         color = factor(word, levels = unique(word))
#       )
#     ) +
#     ggplot2::scale_size_continuous(range = c(3, 8)) +
#     ggplot2::scale_color_manual(values = colors) +
#     ggplot2::theme_void() +
#     ggplot2::theme(
#       legend.position = "none",
#       panel.background = ggplot2::element_rect(fill = bg_color),
#       ...
#     )
  
#   return(p)
# }

# # # Example usage:
# # sentences <- c("This is a test sentence.", "Another sentence for testing.")
# # wordcloud_generator(sentences)
