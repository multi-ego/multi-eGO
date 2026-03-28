import { Routes, Route } from "react-router-dom";
import Navbar from "./components/Navbar.jsx";
import Footer from "./components/Footer.jsx";
import Landing from "./pages/Landing.jsx";
import ConfigBuilder from "./pages/ConfigBuilder.jsx";
import SetupAssistant from "./pages/SetupAssistant.jsx";

export default function App() {
  return (
    <div className="flex min-h-screen flex-col">
      <Navbar />
      <main className="flex-1">
        <Routes>
          <Route path="/" element={<Landing />} />
          <Route path="/setup" element={<SetupAssistant />} />
          <Route path="/config" element={<ConfigBuilder />} />
        </Routes>
      </main>
      <Footer />
    </div>
  );
}
