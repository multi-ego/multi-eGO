import { defineConfig } from "vite";
import react from "@vitejs/plugin-react";

export default defineConfig({
  plugins: [react()],
  // Set base to repo name for GitHub Pages — update if the repo is renamed
  base: "/multi-eGO/",
});
