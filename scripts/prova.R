pdf("plot.pdf")

plot_df <- final %>%
  group_by(block, is_congruent_trial, trials_after_clip) %>%
  summarise(
    y = mean(rtTukey, na.rm = TRUE, trim = 0.05),
    stderr = sqrt(var(rtTukey, na.rm = TRUE) / n())
  )

plot_df$lower <- plot_df$y - plot_df$stderr
plot_df$upper <- plot_df$y + plot_df$stderr

pd <- position_dodge(0.2)

ggplot(plot_df, aes(x = trials_after_clip, y = y, color = is_congruent_trial)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), lwd = 1.05, position = pd) +
  geom_line(position = pd, lwd = 1.05) +
  geom_point(position = pd) +
  ylim(500, 700) +
  facet_wrap(~ block) +
  xlab("Trials after the video clip") +
  ylab("Response Latency (ms)") 
# ggtitle(title)
# geom_hline(yintercept = 0, lty = 2) +
# theme(legend.position = "bottom")

dev.off()
